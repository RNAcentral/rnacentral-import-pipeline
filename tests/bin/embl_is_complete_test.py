import subprocess
from pathlib import Path

BIN = Path(__file__).resolve().parents[2] / "bin" / "embl-is-complete"
EMBL = Path("data/ensembl/Plasmodium_falciparum_3D7.ncrna.embl")


def run(path):
    return subprocess.run([str(BIN), str(path)], capture_output=True, text=True)


def test_accepts_a_whole_embl_file(tmp_path):
    path = tmp_path / "whole.dat"
    path.write_text(EMBL.read_text())
    assert run(path).returncode == 0


# Oxyria digyna's geneset on the Ensembl FTP stops mid-sequence, with no closing '//'.
def test_rejects_an_embl_file_cut_off_mid_record(tmp_path):
    text = EMBL.read_text()
    path = tmp_path / "cut.dat"
    path.write_text(text[: text.rindex("//")])
    result = run(path)
    assert result.returncode != 0
    assert "truncated" in result.stderr
