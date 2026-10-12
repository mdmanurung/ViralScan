# M-01: how `kb count -x 10XV3` strands reads, and what that means for AUD-01 grouping

Read-only investigation, 2026-10-11 (bundle P2). Nothing was run against data; no `kb ref` / `kb count`.

## What was checked

| Source | Finding |
|---|---|
| `kb_python/count.py:220-225` (kb-python 0.30.2) | `kallisto bus` receives `--unstranded`, `--fr-stranded` or `--rf-stranded` **only** when `strand` is set; with `strand=None` no strand flag is passed. |
| `kb_python/main.py:436-452` | A strand is inferred only from a `-x` string with a `%` suffix (`-x 10XV3%FORWARD`). Plain `-x 10XV3` leaves `args.strand = None`. |
| `kb --list` | `10XV3` is a plain on-list technology (barcode `0,0,16`, UMI `0,16,28`, cDNA read 1); the table carries no strand column. |
| `kallisto bus` help (both the conda 0.52.0 binary and kb's bundled 0.52.0, bustools 0.45.1) | Options `--unstranded`, `--fr-stranded` ("first read forward"), `--rf-stranded` ("first read reverse"). The help does **not** state the default when none is given. |
| `src/viralscan/Snakefile:150,170-176` | `--strand` is passed to `kb count` only when `config["strand"]` is set (`${ST:+--strand "$ST"}`); default `None` omits it. |
| `src/viralscan/menu.py:1113-1126,1318-1320` + `src/viralscan/strand.py` | `--strand auto` runs three 1M-pair pilots (`forward`, `reverse`, `unstranded`), picks `reverse`/`forward` if it keeps >= 0.8 of the unstranded pseudoalignment rate, else `unstranded`; the choice then replaces `args.strand`. Opt-in; recorded as `strand_inference`. |
| `src/viralscan/cli/parsers.py:833-848`, `docs/cli_reference.md:574`, PLAN `DEF-02`/F-020 | Repo's own documented conclusion: with the pinned kallisto 0.52.0, no flag means effectively **unstranded** for every chemistry (10x technologies are marked strand-specific but the strand is left unset; kallisto <= 0.51.1 defaulted 10x to forward). |

## Conclusion

* Default viralscan runs (no `--strand`, `-x 10XV3`) are **unstranded** under kallisto 0.52.0: `kb count` passes no strand flag and the pinned kallisto then treats reads as non-strand-specific. This rests on the repo's F-020 / DEF-02 measurements (HHV-6B, `-x 10xv2`, no flag, kallisto 0.52.0: 60.2 % pseudoaligned, i.e. unstranded behaviour; covid `-x 10xv3`, kallisto **0.51.1**, no flag: 6.4 / 8.9 %, i.e. forward). **Unverified by me:** the kallisto 0.52.0 source default (not in the help text) and a run on this index; I did not execute kallisto on data.
* Explicit `--strand forward|reverse` runs are stranded (`--fr-stranded` / `--rf-stranded`); `auto` picks per sample.
* **Antisense overlaps are therefore ambiguous for the default libraries**: the index holds only the sense sequence of each transcript, and an unstranded read from the antisense strand of gene X that matches gene Y's sense sequence pseudoaligns to Y as well (equivalence class {X, Y}); `bustools count` then drops multi-gene classes from the unique layer unless `--mm`. In an explicitly stranded run, only same-orientation matches would be ambiguous, so the unstranded definition is a superset (conservative).

## Decision for AUD-01

Group **regardless of strand**: genes are compared by *canonical* 31-mers (min of forward and reverse complement), which is the unstranded ambiguity relation and the conservative bound for stranded runs. Grouping is keyed by sequence (and, in the coordinate fallback, by `seqname`), not by strand. A strand-aware variant would be needed only if the catalogue were ever to be valid for `--strand forward|reverse` runs alone; that would shrink some groups (antisense-only pairs would split), never grow them.

Open items (not done here): confirm the kallisto 0.52.0 no-flag default from source or a tiny `kallisto bus` run (needs an index; excluded from this bundle), and whether `docs/panel_integrity.md` / the quickstart should say "unstranded by default" as the CLI help already does.
