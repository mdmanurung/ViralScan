"""Original upstream (mahfouzlab/ViralScan) invocation forms parse like the fork-native forms.

The original used ``type=bool`` value-style flags (``-v False``, ``-ref True``); the fork keeps
bare flags and ``--no-*`` and additionally accepts a strictly validated optional value.
"""

import pytest

from viralscan.menu import _build_run_config, build_parser

README_EXISTING = [
    "-t",
    "transcripts.txt",
    "-i",
    "index.idx",
    "-o",
    "output/",
    "-s1",
    "a_1.fastq.gz",
    "-s2",
    "a_2.fastq.gz",
]
README_BUILD = [
    "-reference",
    "True",
    "-fasta",
    "g.fa",
    "-gtf",
    "g.gtf",
    "-o",
    "output/",
    "-s1",
    "a_1.fastq.gz",
    "-s2",
    "a_2.fastq.gz",
]
KEYS = (
    "output",
    "sample1",
    "sample2",
    "transcripts",
    "index",
    "cores",
    "reference",
    "gtf",
    "fasta",
    "f1",
    "visual",
    "multimapping",
    "umap",
    "technology",
    "whitelist",
)


def _ns(argv):
    return vars(build_parser().parse_args(argv))


def _same(orig, native):
    a, b = _ns(orig), _ns(native)
    for k in KEYS:
        assert a[k] == b[k], k


def test_readme_existing_index_defaults():
    ns = _ns(README_EXISTING)
    assert (ns["transcripts"], ns["index"], ns["reference"], ns["visual"], ns["multimapping"]) == (
        "transcripts.txt",
        "index.idx",
        False,
        True,
        True,
    )


def test_readme_build_index_single_dash_reference():
    _same(README_BUILD, ["--reference"] + README_BUILD[2:])
    _same(README_BUILD, ["-ref"] + README_BUILD[2:])
    assert _ns(README_BUILD)["reference"] is True


def test_reference_bare_flag_before_option_does_not_swallow_it():
    ns = _ns(["-ref", "-fasta", "g.fa", "-gtf", "g.gtf"])
    assert ns["reference"] is True and ns["fasta"] == "g.fa"


@pytest.mark.parametrize(
    "flag,dest,off",
    [
        ("-v", "visual", "--no-visual"),
        ("--visual", "visual", "--no-visual"),
        ("-mm", "multimapping", "--no-multimapping"),
        ("--multimapping", "multimapping", "--no-multimapping"),
    ],
)
@pytest.mark.parametrize(
    "word,val",
    [
        ("False", False),
        ("false", False),
        ("0", False),
        ("no", False),
        ("True", True),
        ("true", True),
        ("1", True),
        ("yes", True),
    ],
)
def test_value_style_booleans(flag, dest, off, word, val):
    # 'False' must be False (the type=bool bug made every non-empty string truthy)
    assert _ns(README_EXISTING + [flag, word])[dest] is val
    assert _ns(README_EXISTING + [flag])[dest] is True
    assert _ns(README_EXISTING + [off])[dest] is False


def test_reference_false_is_false():
    assert _ns(["-reference", "False"])["reference"] is False


@pytest.mark.parametrize(
    "argv", [["-v", "maybe"], ["-mm", "2"], ["-ref", "Flase"], ["--no-visual", "True"]]
)
def test_invalid_value_is_argparse_error(argv):
    with pytest.raises(SystemExit) as e:
        build_parser().parse_args(argv)
    assert e.value.code == 2


@pytest.mark.parametrize("argv", [["-v", "data"], ["-ref", "build-ref"], ["-mm", "data"]])
def test_bare_flag_before_subcommand_takes_it_as_value(argv):
    # Pins documented behaviour: use the flag after the subcommand or `--visual=false`.
    with pytest.raises(SystemExit) as e:
        build_parser().parse_args(argv)
    assert e.value.code == 2


def test_comma_lists_f1_and_misc_flags():
    ns = _ns(
        [
            "-o",
            "out",
            "-s1",
            "a_1.fq,b_1.fq",
            "-s2",
            "a_2.fq,b_2.fq",
            "-t",
            "t.txt",
            "-i",
            "i.idx",
            "-c",
            "4",
            "-f1",
            "c.fa",
            "-x",
            "10xv3",
            "-w",
            "wl.txt",
            "-umap",
        ]
    )
    assert ns["sample1"] == "a_1.fq,b_1.fq" and ns["sample2"] == "a_2.fq,b_2.fq"
    assert (ns["cores"], ns["f1"], ns["technology"], ns["whitelist"], ns["umap"]) == (
        4,
        "c.fa",
        "10xv3",
        "wl.txt",
        True,
    )


def test_original_and_native_give_same_runconfig():
    def cfg(argv):
        a = build_parser().parse_args(argv)
        return _build_run_config(a, "out", "i.idx", "t.txt", None, "a_1.fq", "a_2.fq")

    orig = cfg(README_EXISTING + ["-v", "False", "-mm", "False"])
    native = cfg(README_EXISTING + ["--no-visual", "--no-multimapping"])
    assert orig == native
    assert orig.visual is False and orig.multimapping is False
