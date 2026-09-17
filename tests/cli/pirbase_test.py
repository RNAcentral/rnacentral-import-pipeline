from pathlib import Path

from click.testing import CliRunner

from rnacentral_pipeline.cli import pirbase as cli
from rnacentral_pipeline.databases.pirbase import parser


def test_full_parse_with_nothing_known_writes_empty_output(tmp_path):
    md5s = tmp_path / "md5s"
    md5s.write_text("")
    known = tmp_path / "known.sqlite"
    parser.build_known(md5s, known)
    out = tmp_path / "out"
    out.mkdir()

    result = CliRunner().invoke(
        cli.cli,
        ["parse", "--known", str(known), "cel", "data/pirbase/cel.fa", str(out)],
    )

    assert result.exit_code == 0, result.output
