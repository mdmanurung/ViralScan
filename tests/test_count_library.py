import importlib.util
from pathlib import Path

import pytest

SCRIPT = Path(__file__).resolve().parents[1] / "analysis/translated_search/count_library.py"
spec = importlib.util.spec_from_file_location("count_library", SCRIPT)
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)

CFG = {
    "kb": "kb",
    "kallisto": "kallisto",
    "technology": "10xv3",
    "samples": {"s1": ["r1.fq.gz", "r2.fq.gz"]},
}


def test_command_matches_snakefile_count_rule():
    cmd = m.command(CFG, "s1", Path("i.idx"), Path("g.t2g"), Path("out/s1"), 4, True)
    assert cmd[:7] == ["kb", "count", "--kallisto", "kallisto", "--aa", "--parity", "single"]
    assert "--num" in cmd and "-w" not in cmd and cmd[-2:] == ["r1.fq.gz", "r2.fq.gz"]
    assert cmd[cmd.index("-x") + 1] == "10xv3" and cmd[cmd.index("-t") + 1] == "4"
    assert cmd[cmd.index("--tmp") + 1] == str(Path("out/s1_tmp"))


def test_whitelist_and_no_num():
    cmd = m.command({**CFG, "whitelist": "wl.txt"}, "s1", Path("i"), Path("g"), Path("o"), 2, False)
    assert cmd[cmd.index("-w") + 1] == "wl.txt" and "--num" not in cmd
    with pytest.raises(SystemExit):
        m.command(CFG, "nope", Path("i"), Path("g"), Path("o"), 2, False)
