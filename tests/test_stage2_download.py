import json
from datetime import timedelta

import pytest

from stage2_helpers import MockArchive, archive_key, make_tart_hdf5, utc
from tart_gnss.stage2.archive import find_observations
from tart_gnss.stage2.download import download_observation
from tart_gnss.stage2.errors import DownloadError
from tart_gnss.stations import Station
from tart_gnss.tle import sha256_file

ALPHA = Station("alpha", "a", -20.0, 57.0, 0.0, "test", tart_id="aa-alpha", archive_name="aa-alpha")
T = utc(2026, 6, 25, 10, 1)


@pytest.fixture
def setup(tmp_path):
    arch = MockArchive()
    arch.add_file(archive_key("aa-alpha", T), make_tart_hdf5(tmp_path / "src.hdf", "alpha", T, n=3))
    rec = find_observations(ALPHA, T - timedelta(minutes=1), T, arch)[0]
    return arch, rec, tmp_path / "raw"


def test_sha256_deterministic(tmp_path):
    import hashlib

    p = tmp_path / "f.bin"
    p.write_bytes(b"tart")
    assert sha256_file(p) == sha256_file(p) == hashlib.sha256(b"tart").hexdigest()


def test_download_records_provenance(setup):
    arch, rec, raw = setup
    f = download_observation(rec, raw, arch)
    assert f.original_filename == rec.filename and f.md5_matches_etag is True
    assert f.sha256 == sha256_file(f.local_path) and f.size_bytes == rec.size_bytes
    prov = json.loads((raw / (rec.filename + ".prov.json")).read_text())
    assert prov["remote_identifier"] == rec.remote_identifier and prov["retrieval_utc"].endswith("Z")


def test_verified_reuse_without_download(setup):
    arch, rec, raw = setup
    download_observation(rec, raw, arch)
    again = download_observation(rec, raw, arch)
    assert again.reused and len(arch.downloads) == 1


def test_tampered_file_refused(setup):
    arch, rec, raw = setup
    f = download_observation(rec, raw, arch)
    with open(f.local_path, "ab") as h:
        h.write(b"x")
    with pytest.raises(DownloadError, match="SHA-256"):
        download_observation(rec, raw, arch)


def test_unknown_existing_file_never_overwritten(setup):
    arch, rec, raw = setup
    raw.mkdir(parents=True)
    (raw / rec.filename).write_bytes(b"something else")
    with pytest.raises(DownloadError, match="without provenance"):
        download_observation(rec, raw, arch)
    assert (raw / rec.filename).read_bytes() == b"something else"


def test_corrupted_transfer_detected(setup):
    arch, rec, raw = setup
    data = bytearray(arch.objects[rec.remote_identifier])
    data[len(data) // 2] ^= 0xFF  # flip one byte, same size
    arch.objects[rec.remote_identifier] = bytes(data)
    with pytest.raises(DownloadError, match="ETag"):
        download_observation(rec, raw, arch)
    assert not (raw / rec.filename).exists()
