# M-01 audit: strand handling across count, replay and read extraction

Date 2026-10-12, branch `codex/aud-12a`, PLAN `AUD-12`. Companion to `M01_strand_note.md` (how `kb count` strands reads, for AUD-01 grouping); this note covers whether every path that re-reads the run uses the strand the count used. Biomni's wording was "kb count defaults 10xv3 forward-stranded but kb extract defaults differ".

## Short answer

No path silently uses a different strand from the count that produced the run. The one real count-vs-replay mismatch was the `-x` technology (not the strand): fixed. Three smaller items are recorded as design or remaining.

## What was checked

| Question | Evidence | Result |
|---|---|---|
| Does ViralScan call `kb extract`? | `grep` of `src/`, `scripts/`, `docs/` for `kb extract` / `"extract"` | **No.** Read extraction is `viralscan evidence` -> `evidence_run.run_evidence`, which re-runs `kallisto bus -n` plus `bustools capture` itself. `kb extract --strand` (help: "default: 'unstranded'", always `technology="bulk"`, single FASTQ) is never reached, so the premise of M-01 does not apply to this code. |
| What does count pass? | `Snakefile:150,182-188` | `--strand` only when `config.strand` is set; kb-python 0.30.2 `count.py:220-225` turns it into `--unstranded/--fr-stranded/--rf-stranded`; unset adds no flag. |
| What does the replay pass? | `evidence.py` `_KALLISTO_STRAND_FLAGS`, `replay_exact_target_bus`; `evidence_run.py:396-412` | Same mapping from `config.strand`; unset adds no flag. Equal to count for all four values. |
| Where is the strand recorded? | `menu.py:1329-1332`, `_build_run_config`, `RunConfig.strand`; `run_safety.py` | `--strand auto` is resolved (`_resolve_auto_strand`) into `sample_args` before `_build_run_config`, so the per-sample `config.yaml` always holds a concrete strand or null, never `auto`. `run_manifest.json` keeps `options.strand` (`auto` or explicit) and, for `auto`, `strand_inference[sample]` with rates and choice (outside the fingerprint). The evidence replay reads `config.yaml`, so it uses the resolved strand. A hand-edited `auto` falls back to no flag with a warning. |
| Do `rerun-multimap` / `rerun-programs` / `auto_evidence` change it? | `menu.py:242-395`, `scripts/auto_evidence.py`, `reads.py` | They never run kallisto on reads except through `run_evidence` (config strand). `rerun-multimap` copies `config.yaml` and keeps its mtime, so `kb_count` is not re-run and the strand is unchanged. `multimap.py` only runs `bustools` on the existing BUS. `extract_virus_reads` builds a Namespace for `run_evidence`; it has no strand parameter and cannot override the config. |
| Is "no flag" the same thing in count and replay? | Live check, kallisto 0.52.0 on a synthetic 10x-style library (half sense, half antisense reads), `tests/integration/test_strand_replay_kallisto.py` | `-x 10xv3` with no flag pseudoaligns 200/200 = `--unstranded`; `--fr-stranded` keeps 100/200. So under the pinned kallisto an unset strand is unstranded on both sides. This closes the open item in `M01_strand_note.md` (the 0.52.0 no-flag default). It is kallisto-version dependent (0.51.1 defaulted 10x to forward, F-020), which is why the integration test pins it. |
| Same technology string? | kb-python `count.py:184-186`; live kallisto | **Mismatch found.** kb count rewrites `10XV4` to `10XV3` before `kallisto bus`; the replay passed the raw `10xv4`, which kallisto 0.52.0 rejects ("technology string must contain two colons"). `viralscan evidence` on a counted 10xv4 run therefore died. |

## What changed

* `evidence.py`: `kallisto_technology()` applies kb's 10XV4 -> 10XV3 rewrite; `replay_exact_target_bus` uses it for `-x`. Commit `6262e3b`.
* `evidence_run.py`: logs the replay strand ("from config.yaml, what kb count used" or "no flag"), so the effective setting is visible.
* Tests: `tests/test_evidence.py` (technology reaches kallisto as kb count passes it; fails on the old code for 10xv4/10XV4), `tests/test_evidence_subcommand.py` (config strand None/forward/reverse/unstranded reaches the replay; hand-edited `auto` gives no flag), `tests/test_strand_consistency.py` (the strand round-trips `config.yaml`), `tests/integration/test_strand_replay_kallisto.py` (live: no flag == unstranded; raw 10xv4 rejected, rewritten accepted; `pytest -m integration`).

## Design decisions, not changed

* **Unset strand means unstranded under kallisto 0.52.0** (and `kb count` records no strand). That is the documented default (`--strand` help, DEF-02); changing it would change every count. Replay deliberately mirrors "no flag" rather than forcing `--unstranded`, so a future kallisto change moves count and replay together.
* **Antisense ambiguity** for default libraries: see `M01_strand_note.md` (AUD-01 grouping is strand-independent on purpose).

## Remaining

1. `evidence_manifest.json` does not record the replay strand/technology (they are in `config.yaml` and now in the log). Adding a manifest key needs the evidence-manifest schema and `docs/output_reference.md` (pinned) to move together; left for a doc/governance pass.
2. `-x 10xv3%FORWARD` (kb's strand-suffix syntax, reachable only with `--force-technology`): kb count strands the run, but `config.strand` stays null, and the replay would hand kallisto the raw `%` string and fail. Not changed; the clean fix is to reject `%` in `-x` or split it into `--strand` in `menu.py`, a behaviour choice for the maintainer.
3. `--strand auto` pilot (`strand.py::_kb_pilot`) runs on the raw FASTQs and passes the raw `-w` (a `.gz` list is decompressed for the real count in the Snakefile but not for the pilot; `p_pseudoaligned` is computed before barcode correction, so `-w` should not change the rate, but a gzipped list may make the pilot fail). Strandedness is a property of the library, so raw vs host/read-filtered reads is not an inconsistency. I could not run `kb count` end to end here (kb's bundled kallisto hung on its input pipes in this sandbox), so the gz-`-w` pilot failure is unreproduced.
4. `chemistry_check.sanity_gate` message calls the low host-index rate "the forward-strand default (F-020)"; true for kallisto <= 0.51.1 only. Wording only.
