import hashlib
import os

import pytest
import requests

from funannotate import test as fun_test


def _fake_download(content_by_url, calls):
    def fake(url, name):
        calls.append(url)
        content = content_by_url.get(url)
        if content is None:
            raise requests.exceptions.RequestException("no such url")
        with open(name, "wb") as f:
            f.write(content)

    return fake


def test_every_key_has_release_name_osf_link_and_checksum():
    assert set(fun_test.TEST_DATA_LINKS) == set(fun_test.TEST_DATA_OSF)
    for name in fun_test.TEST_DATA_LINKS.values():
        assert len(fun_test.TEST_DATA_SHA256[name]) == 64


def test_download_from_release_verifies_checksum(tmp_path, monkeypatch):
    data = b"test archive"
    monkeypatch.setitem(fun_test.TEST_DATA_SHA256, "test-mask.tar.gz", hashlib.sha256(data).hexdigest())
    calls = []
    url = fun_test.TEST_DATA_URL.format("test-mask.tar.gz")
    monkeypatch.setattr(fun_test, "download", _fake_download({url: data}, calls))
    monkeypatch.chdir(tmp_path)
    fun_test.download_test_data("mask", "out.tar.gz")
    assert calls == [url]
    assert (tmp_path / "out.tar.gz").read_bytes() == data


def test_falls_back_to_osf_when_release_fails(tmp_path, monkeypatch):
    data = b"osf copy"
    monkeypatch.setitem(fun_test.TEST_DATA_SHA256, "test-mask.tar.gz", hashlib.sha256(data).hexdigest())
    calls = []
    osf = fun_test.TEST_DATA_OSF["mask"]
    monkeypatch.setattr(fun_test, "download", _fake_download({osf: data}, calls))
    monkeypatch.chdir(tmp_path)
    fun_test.download_test_data("mask", "out.tar.gz")
    assert calls == [fun_test.TEST_DATA_URL.format("test-mask.tar.gz"), osf]


def test_checksum_mismatch_everywhere_exits_and_removes_file(tmp_path, monkeypatch):
    calls = []
    urls = {fun_test.TEST_DATA_URL.format("test-mask.tar.gz"): b"bad", fun_test.TEST_DATA_OSF["mask"]: b"bad"}
    monkeypatch.setattr(fun_test, "download", _fake_download(urls, calls))
    monkeypatch.chdir(tmp_path)
    with pytest.raises(SystemExit):
        fun_test.download_test_data("mask", "out.tar.gz")
    assert not os.path.exists("out.tar.gz")
