#!/usr/bin/env python3
"""Mirror funannotate's BUSCO lineage tarballs from osf.io to a GitHub release.

`funannotate setup -b ...` downloads the 28 BUSCO (odb9) lineage tarballs listed
under "busco" in funannotate/downloads.json. They are hosted on the Funannotate
OSF project (osf.io/bj7v4), which serves them at roughly 0.5 MB/s, so a full
`-b all` takes ~2 h for 3.3 GB. This script copies them verbatim to assets on a
GitHub release of this repo and rewrites downloads.json to point there.

setupDB.py fetches downloads.json from master at runtime, so merging the
rewritten JSON switches every install over without a new funannotate release.

Steps (each is idempotent; re-run to resume):
  download   fetch each tarball from osf.io into --dir, verify it against the
             sha256 OSF publishes, and check its top-level folder is the name
             downloads.json expects (setupDB renames that folder after untar)
  release    create the GitHub release (--tag) if missing and upload any asset
             not already there (needs `gh` authenticated with push access)
  json       rewrite the "busco" URLs in downloads.json to the release assets
  verify     HEAD every URL in downloads.json "busco" and compare its size with
             the OSF size (run after `json`, once the release is published)

Usage:
  python scripts/mirror_busco_to_github.py download --dir /tmp/busco_mirror
  python scripts/mirror_busco_to_github.py release  --dir /tmp/busco_mirror
  python scripts/mirror_busco_to_github.py json
  python scripts/mirror_busco_to_github.py verify
"""

import argparse
import concurrent.futures
import hashlib
import json
import os
import re
import subprocess
import sys
import tarfile

import requests

REPO = "nextgenusfs/funannotate"
DEFAULT_TAG = "busco-odb9-mirror"
HERE = os.path.dirname(os.path.abspath(__file__))
DOWNLOADS_JSON = os.path.join(HERE, "..", "funannotate", "downloads.json")

RELEASE_NOTES = """\
Verbatim mirror of the 28 BUSCO v2 (OrthoDB v9, "odb9") lineage datasets that
`funannotate setup -b ...` installs, previously served only from the Funannotate
OSF project (https://osf.io/bj7v4/). Hosted here for faster, more reliable
downloads; the files are byte-identical (sha256 verified against OSF, listed in
SHA256SUMS).

Dataset credit: BUSCO lineage datasets, Simão et al. 2015 (Bioinformatics
31:3210) and Waterhouse et al. 2018 (Mol Biol Evol 35:543), built from OrthoDB
v9 (Zdobnov et al. 2017, Nucleic Acids Res 45:D744). https://busco.ezlab.org/

This is a data-only release, not a funannotate software release.
"""


def load_busco():
    with open(DOWNLOADS_JSON) as fh:
        return json.load(fh)["busco"]


def osf_id(url):
    m = re.match(r"https://osf\.io/([a-z0-9]+)/", url)
    return m.group(1) if m else None


def osf_meta(fid):
    r = requests.get("https://api.osf.io/v2/files/%s/" % fid, timeout=60)
    r.raise_for_status()
    a = r.json()["data"]["attributes"]
    return {"name": a["name"], "size": a["size"], "sha256": a["extra"]["hashes"]["sha256"]}


def sha256sum(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def top_level_dirs(path):
    with tarfile.open(path, "r:gz") as tf:
        return {m.name.lstrip("./").split("/")[0] for m in tf.getmembers() if m.name.strip("./")}


def fetch_one(key, url, folder, outdir):
    fid = osf_id(url)
    if not fid:
        return key, None, "not an osf.io URL (already mirrored?): %s" % url
    meta = osf_meta(fid)
    dest = os.path.join(outdir, meta["name"])
    if not (os.path.isfile(dest) and os.path.getsize(dest) == meta["size"] and sha256sum(dest) == meta["sha256"]):
        tmp = dest + ".part"
        with requests.get(url, stream=True, allow_redirects=True, timeout=120) as r:
            r.raise_for_status()
            with open(tmp, "wb") as fh:
                for chunk in r.iter_content(chunk_size=1 << 20):
                    fh.write(chunk)
        if sha256sum(tmp) != meta["sha256"]:
            os.remove(tmp)
            return key, None, "sha256 mismatch vs OSF"
        os.replace(tmp, dest)
    tops = top_level_dirs(dest)
    if tops != {folder}:
        return key, None, "tarball top-level %s != expected %s" % (sorted(tops), folder)
    meta.update(key=key, folder=folder, osf_url=url)
    return key, meta, None


def cmd_download(args):
    os.makedirs(args.dir, exist_ok=True)
    busco = load_busco()
    results, errors = {}, {}
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.jobs) as ex:
        futs = [ex.submit(fetch_one, k, v[0], v[1], args.dir) for k, v in busco.items()]
        for f in concurrent.futures.as_completed(futs):
            key, meta, err = f.result()
            if err:
                errors[key] = err
                print("FAIL %-26s %s" % (key, err), flush=True)
            else:
                results[key] = meta
                print("ok   %-26s %-32s %8.1f MB" % (key, meta["name"], meta["size"] / 1e6), flush=True)
    with open(os.path.join(args.dir, "manifest.json"), "w") as fh:
        json.dump(results, fh, indent=2, sort_keys=True)
    with open(os.path.join(args.dir, "SHA256SUMS"), "w") as fh:
        for m in sorted(results.values(), key=lambda m: m["name"]):
            fh.write("%s  %s\n" % (m["sha256"], m["name"]))
    print("%d ok, %d failed, %.2f GB" % (len(results), len(errors), sum(m["size"] for m in results.values()) / 1e9))
    return 1 if errors else 0


def gh(*a, check=True):
    return subprocess.run(["gh", *a], check=check, capture_output=True, text=True)


def cmd_release(args):
    with open(os.path.join(args.dir, "manifest.json")) as fh:
        manifest = json.load(fh)
    if len(manifest) != len(load_busco()):
        sys.exit("manifest has %d entries, downloads.json has %d; re-run download" % (len(manifest), len(load_busco())))
    if gh("release", "view", args.tag, "-R", REPO, check=False).returncode != 0:
        gh("release", "create", args.tag, "-R", REPO, "--title", "BUSCO odb9 lineage mirror",
           "--notes", RELEASE_NOTES, "--latest=false")
        print("created release %s" % args.tag)
    assets = json.loads(gh("release", "view", args.tag, "-R", REPO, "--json", "assets").stdout)["assets"]
    have = {a["name"] for a in assets}
    files = [m["name"] for m in manifest.values()] + ["SHA256SUMS"]
    todo = [f for f in files if f not in have]
    for f in todo:
        gh("release", "upload", args.tag, "-R", REPO, os.path.join(args.dir, f))
        print("uploaded %s" % f, flush=True)
    print("%d uploaded, %d already present" % (len(todo), len(files) - len(todo)))


def cmd_json(args):
    with open(os.path.join(args.dir, "manifest.json")) as fh:
        manifest = json.load(fh)
    with open(DOWNLOADS_JSON) as fh:
        data = json.load(fh)
    base = "https://github.com/%s/releases/download/%s/" % (REPO, args.tag)
    for key, (url, folder) in data["busco"].items():
        data["busco"][key] = [base + manifest[key]["name"], folder]
    with open(DOWNLOADS_JSON, "w") as fh:
        json.dump(data, fh, indent=2)
        fh.write("\n")
    print("rewrote %d busco URLs -> %s" % (len(data["busco"]), base))


def cmd_verify(args):
    manifest = {}
    mpath = os.path.join(args.dir, "manifest.json")
    if os.path.isfile(mpath):
        with open(mpath) as fh:
            manifest = json.load(fh)
    bad = 0
    for key, (url, folder) in load_busco().items():
        r = requests.head(url, allow_redirects=True, timeout=60)
        size = int(r.headers.get("content-length", -1))
        want = manifest.get(key, {}).get("size")
        ok = r.status_code == 200 and (want is None or size == want)
        bad += not ok
        print("%s %-26s %s %d %s" % ("ok  " if ok else "FAIL", key, r.status_code, size, url))
    print("%d/%d reachable with expected size" % (len(load_busco()) - bad, len(load_busco())))
    return 1 if bad else 0


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("step", choices=["download", "release", "json", "verify"])
    p.add_argument("--dir", default="busco_mirror", help="local staging directory")
    p.add_argument("--tag", default=DEFAULT_TAG, help="GitHub release tag")
    p.add_argument("--jobs", type=int, default=4, help="parallel downloads")
    args = p.parse_args()
    return {"download": cmd_download, "release": cmd_release, "json": cmd_json, "verify": cmd_verify}[args.step](args) or 0


if __name__ == "__main__":
    sys.exit(main())
