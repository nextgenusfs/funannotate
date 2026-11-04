"""The install scripts pin Trinity, PASA and EVM to release tags.

A pin that is a branch name moves whenever that branch is pushed to, so a
rebuilt container or conda environment can silently change. Each script takes a
default ref from  ${NAME_RUST_COMMIT:-<ref>}; that ref must be a tag like v2.5.0.
"""
import os
import re

import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PINS = {
    "pixi_install_evm.sh": "EVM_RUST_COMMIT",
    "pixi_install_pasa.sh": "PASA_RUST_COMMIT",
    "pixi_install_trinity.sh": "TRINITY_RUST_COMMIT",
}
TAG = re.compile(r"^v\d+\.\d+\.\d+([._-][A-Za-z0-9.]+)?$")


def default_ref(script, variable):
    text = open(os.path.join(ROOT, "install_scripts", script)).read()
    match = re.search(r"\$\{%s:-([^}]+)\}" % variable, text)
    assert match, "%s has no default for %s" % (script, variable)
    return match.group(1)


@pytest.mark.parametrize("script,variable", sorted(PINS.items()))
def test_default_ref_is_a_release_tag(script, variable):
    ref = default_ref(script, variable)
    assert TAG.match(ref), "%s defaults to %r, which is not a vN.N.N release tag" % (script, ref)


def version_of(ref):
    return tuple(int(x) for x in re.match(r"v(\d+)\.(\d+)\.(\d+)", ref).groups())


def test_evm_pin_has_the_options_runevm_uses():
    """funannotate-runEVM.py passes --threads, --skip_empty and --concat_out when the Rust EVM tools list
    them (evm_cmds.py); they first appeared in v2.5.1_rust, so an older pin would leave them unused."""
    assert version_of(default_ref("pixi_install_evm.sh", "EVM_RUST_COMMIT")) >= (2, 5, 1)

