#!/usr/bin/env python3
"""Mirror funannotate's static database tarballs to a GitHub release.

`funannotate setup` downloads several static tarballs that were only hosted on
the Funannotate OSF project (osf.io/bj7v4):
  - the 28 BUSCO (odb9) lineage tarballs under "busco" in downloads.json
  - the repeat protein library ("repeats") and the pre-computed BUSCO
    outgroups ("outgroups") under "downloads"
OSF serves them slowly (~0.5 MB/s) and through a 308 redirect that urllib
before Python 3.11 cannot follow. This script copies them verbatim to assets on
a GitHub release of this repo and rewrites downloads.json to point there.

Files must stay byte-identical: `setup -u` compares the remote md5 with the md5
stored in each user's database, so a changed file would force a re-download.

setupDB.py fetches downloads.json from master at runtime, so merging the
rewritten JSON switches every install over without a new funannotate release.

A source URL may be an osf.io link (verified against the sha256 OSF publishes)
or an asset of an earlier GitHub mirror release (verified against that
release's SHA256SUMS), so the script can also move files between releases.

Steps (each is idempotent; re-run to resume):
  download   fetch each tarball into --dir, verify its sha256, and check its
             top-level entry is the name setupDB expects after untar
  release    create the GitHub release (--tag) if missing and upload any asset
             not already there (needs `gh` authenticated with push access)
  json       rewrite the mirrored URLs in downloads.json to the release assets
  verify     HEAD every mirrored URL in downloads.json and compare its size
             with the manifest (run after `json`, once the release is published)

Usage:
  python scripts/mirror_db_to_github.py download --dir /tmp/db_mirror
  python scripts/mirror_db_to_github.py release  --dir /tmp/db_mirror
  python scripts/mirror_db_to_github.py json     --dir /tmp/db_mirror
  python scripts/mirror_db_to_github.py verify   --dir /tmp/db_mirror
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
DEFAULT_TAG = "db-mirror"
HERE = os.path.dirname(os.path.abspath(__file__))
DOWNLOADS_JSON = os.path.join(HERE, "..", "funannotate", "downloads.json")

# "downloads" entries to mirror -> top-level tar entry setupDB.py expects
STATIC_DOWNLOADS = {
    "repeats": "funannotate.repeat.proteins.fa",
    "outgroups": "outgroups",
}

RELEASE_TITLE = "funannotate database mirror"
RELEASE_NOTES = """\
Verbatim mirror of the static database tarballs that `funannotate setup`
installs, previously served only from the Funannotate OSF project
(https://osf.io/bj7v4/). Hosted here for faster, more reliable downloads; the
files are byte-identical (sha256 verified against OSF, listed in SHA256SUMS).

Contents:
- 28 BUSCO v2 (OrthoDB v9, "odb9") lineage datasets (`funannotate setup -b`)
- funannotate.repeat.proteins.fa.tar.gz: repeat protein library used to build
  the repeats DIAMOND database
- busco_outgroups.tar.gz: pre-computed BUSCO (dikarya) protein sets for six
  outgroup species used by `funannotate compare`

BUSCO dataset credit: Simão et al. 2015 (Bioinformatics 31:3210) and
Waterhouse et al. 2018 (Mol Biol Evol 35:543), built from OrthoDB v9
(Zdobnov et al. 2017, Nucleic Acids Res 45:D744). https://busco.ezlab.org/

This is a data-only release, not a funannotate software release.
"""


def load_json():
    with open(DOWNLOADS_JSON) as fh:
        return json.load(fh)


def load_entries():
    """Return {manifest key: (url, expected top-level tar entry)}."""
    data = load_json()
    entries = {"busco/%s" % k: (v[0], v[1]) for k, v in data["busco"].items()}
    for k, top in STATIC_DOWNLOADS.items():
        entries["downloads/%s" % k] = (data["downloads"][k], top)
    return entries


def osf_id(url):
    m = re.match(r"https://osf\.io/([a-z0-9]+)/", url)
    return m.group(1) if m else None


def gh_asset(url):
    m = re.match(r"https://github\.com/([^/]+/[^/]+)/releases/download/([^/]+)/([^/]+)$", url)
    return m.groups() if m else None


def osf_meta(fid):
    r = requests.get("https://api.osf.io/v2/files/%s/" % fid, timeout=60)
    r.raise_for_status()
    a = r.json()["data"]["attributes"]
    return {"name": a["name"], "sha256": a["extra"]["hashes"]["sha256"]}


_sums_cache = {}


def release_sums(repo, tag):
    if (repo, tag) not in _sums_cache:
        url = "https://github.com/%s/releases/download/%s/SHA256SUMS" % (repo, tag)
        r = requests.get(url, timeout=60)
        r.raise_for_status()
        _sums_cache[(repo, tag)] = {n: s for s, n in (line.split() for line in r.text.splitlines() if line.strip())}
    return _sums_cache[(repo, tag)]


def source_meta(url):
    fid = osf_id(url)
    if fid:
        return osf_meta(fid)
    asset = gh_asset(url)
    if asset:
        repo, tag, name = asset
        sums = release_sums(repo, tag)
        if name not in sums:
            raise ValueError("%s not listed in %s %s SHA256SUMS" % (name, repo, tag))
        return {"name": name, "sha256": sums[name]}
    raise ValueError("unsupported source URL: %s" % url)


def sha256sum(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def md5sum(path):
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def top_level_entries(path):
    with tarfile.open(path, "r:gz") as tf:
        return {m.name.lstrip("./").split("/")[0] for m in tf.getmembers() if m.name.strip("./")}


def fetch_one(key, url, top, outdir):
    try:
        meta = source_meta(url)
    except (ValueError, requests.RequestException) as e:
        return key, None, str(e)
    dest = os.path.join(outdir, meta["name"])
    if not (os.path.isfile(dest) and sha256sum(dest) == meta["sha256"]):
        tmp = dest + ".part"
        with requests.get(url, stream=True, allow_redirects=True, timeout=120) as r:
            r.raise_for_status()
            with open(tmp, "wb") as fh:
                for chunk in r.iter_content(chunk_size=1 << 20):
                    fh.write(chunk)
        if sha256sum(tmp) != meta["sha256"]:
            os.remove(tmp)
            return key, None, "sha256 mismatch vs source"
        os.replace(tmp, dest)
    tops = top_level_entries(dest)
    if tops != {top}:
        return key, None, "tarball top-level %s != expected %s" % (sorted(tops), top)
    meta.update(key=key, top=top, source_url=url, size=os.path.getsize(dest), md5=md5sum(dest))
    return key, meta, None


def cmd_download(args):
    os.makedirs(args.dir, exist_ok=True)
    entries = load_entries()
    results, errors = {}, {}
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.jobs) as ex:
        futs = [ex.submit(fetch_one, k, url, top, args.dir) for k, (url, top) in entries.items()]
        for f in concurrent.futures.as_completed(futs):
            key, meta, err = f.result()
            if err:
                errors[key] = err
                print("FAIL %-36s %s" % (key, err), flush=True)
            else:
                results[key] = meta
                print("ok   %-36s %-40s %8.1f MB" % (key, meta["name"], meta["size"] / 1e6), flush=True)
    with open(os.path.join(args.dir, "manifest.json"), "w") as fh:
        json.dump(results, fh, indent=2, sort_keys=True)
    with open(os.path.join(args.dir, "SHA256SUMS"), "w") as fh:
        for m in sorted(results.values(), key=lambda m: m["name"]):
            fh.write("%s  %s\n" % (m["sha256"], m["name"]))
    print("%d ok, %d failed, %.2f GB" % (len(results), len(errors), sum(m["size"] for m in results.values()) / 1e9))
    return 1 if errors else 0


def gh(*a, check=True):
    return subprocess.run(["gh", *a], check=check, capture_output=True, text=True)


def load_manifest(args, complete=True):
    mpath = os.path.join(args.dir, "manifest.json")
    if not os.path.isfile(mpath):
        if complete:
            sys.exit("no manifest at %s; run download first" % mpath)
        return {}
    with open(mpath) as fh:
        manifest = json.load(fh)
    if complete and set(manifest) != set(load_entries()):
        sys.exit("manifest has %d entries, downloads.json has %d; re-run download" % (len(manifest), len(load_entries())))
    return manifest


def cmd_release(args):
    manifest = load_manifest(args)
    if gh("release", "view", args.tag, "-R", REPO, check=False).returncode != 0:
        gh("release", "create", args.tag, "-R", REPO, "--title", RELEASE_TITLE,
           "--notes", RELEASE_NOTES, "--latest=false")
        print("created release %s" % args.tag)
    assets = json.loads(gh("release", "view", args.tag, "-R", REPO, "--json", "assets").stdout)["assets"]
    have = {a["name"] for a in assets}
    files = [m["name"] for m in manifest.values()]
    todo = [f for f in files if f not in have]
    for f in todo:
        gh("release", "upload", args.tag, "-R", REPO, os.path.join(args.dir, f))
        print("uploaded %s" % f, flush=True)
    # SHA256SUMS changes whenever files are added, so always replace it
    gh("release", "upload", args.tag, "-R", REPO, os.path.join(args.dir, "SHA256SUMS"), "--clobber")
    print("%d uploaded, %d already present, SHA256SUMS refreshed" % (len(todo), len(files) - len(todo)))


def cmd_json(args):
    manifest = load_manifest(args)
    data = load_json()
    base = "https://github.com/%s/releases/download/%s/" % (REPO, args.tag)
    for key, (url, folder) in data["busco"].items():
        data["busco"][key] = [base + manifest["busco/%s" % key]["name"], folder]
    for key in STATIC_DOWNLOADS:
        data["downloads"][key] = base + manifest["downloads/%s" % key]["name"]
    with open(DOWNLOADS_JSON, "w") as fh:
        json.dump(data, fh, indent=2)
        fh.write("\n")
    print("rewrote %d URLs -> %s" % (len(manifest), base))


def cmd_verify(args):
    manifest = load_manifest(args, complete=False)
    entries = load_entries()
    bad = 0
    for key, (url, top) in entries.items():
        r = requests.head(url, allow_redirects=True, timeout=60)
        size = int(r.headers.get("content-length", -1))
        want = manifest.get(key, {}).get("size")
        ok = r.status_code == 200 and (want is None or size == want)
        bad += not ok
        print("%s %-36s %s %d %s" % ("ok  " if ok else "FAIL", key, r.status_code, size, url))
    print("%d/%d reachable with expected size" % (len(entries) - bad, len(entries)))
    return 1 if bad else 0


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("step", choices=["download", "release", "json", "verify"])
    p.add_argument("--dir", default="db_mirror", help="local staging directory")
    p.add_argument("--tag", default=DEFAULT_TAG, help="GitHub release tag")
    p.add_argument("--jobs", type=int, default=4, help="parallel downloads")
    args = p.parse_args()
    return {"download": cmd_download, "release": cmd_release, "json": cmd_json, "verify": cmd_verify}[args.step](args) or 0


if __name__ == "__main__":
    sys.exit(main())
