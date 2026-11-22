import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
CLUSTER = (ROOT / "config" / "cluster.config").read_text()


def test_r2dt_container_tmp_is_node_local():
    # One shared --workdir on /hps put every layout job's /tmp in a single Lustre
    # directory; release 29 ran at 11-16% CPU with traveler stuck in D state.
    options = re.search(r'r2dt_container\s*=\s*"([^"]*)"', CLUSTER).group(1)
    workdir = re.search(r"--workdir\s+(\S+)", options).group(1)
    assert workdir == "/tmp"
