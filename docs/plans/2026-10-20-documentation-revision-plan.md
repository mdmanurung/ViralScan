# Documentation revision plan

Status: implemented on `codex/docs-revision` (2026-10-20, PLAN `DOCS-REV-01`), not yet merged or published. Checks run: strict Sphinx build (`-W --keep-going`, fallback theme because `scanpydoc` is not installed locally), docs consistency/sources tests, ship-scope tests, notebook parse, relative-link check. The `scanpydoc` theme build runs in CI (`docs.yml`).
Reviewed: 2026-10-20. Scope: the published documentation website and its local sources.

## Objective

Give readers one clear path: install ViralScan, run a sample, inspect the results,
and assess the evidence behind a candidate. Keep detailed definitions available
without requiring readers to learn every feature before their first run.

## Findings from the published site

- The home page repeats installation and multiple run commands, while linking
  directly to API material. It has no equally prominent route to assessing results.
- The sidebar lists both `API Reference` and `API reference (autodoc)`. Two
  tutorial links display `[skip-ci]` as their titles.
- The quickstart presents three reference modes before the combined host + virus
  workflow that it recommends for routine host-aware analysis. Optional enrichment,
  D-list construction, and a long flag table distract from finishing the first run.
- The output reference combines file locations, column dictionaries, evidence
  interpretation, optional analyses, implementation history, and study-specific
  observations. The FAQ repeats substantial interpretation and reference guidance.
- The reference-panel page opens with the extended anellovirus reference before
  explaining the general panel. Users need the general choice first.
- The quickstart suggests D-list construction for host homology, while the CLI
  reference contains a study-specific comparison warning against treating it as
  strong host control. Consolidate this guidance, retaining the limits of that
  comparison rather than turning it into a universal performance ranking.

Sources:
[home](https://mdmanurung.github.io/ViralScan/index.html),
[installation](https://mdmanurung.github.io/ViralScan/installation.html),
[quickstart](https://mdmanurung.github.io/ViralScan/quickstart.html),
[outputs](https://mdmanurung.github.io/ViralScan/output_reference.html),
[FAQ](https://mdmanurung.github.io/ViralScan/faq.html),
[CLI](https://mdmanurung.github.io/ViralScan/cli_reference.html),
[panel](https://mdmanurung.github.io/ViralScan/reference_panel.html),
[tutorial index](https://mdmanurung.github.io/ViralScan/vignettes/README.html).
Also inspected the published counting contract, panel integrity, and workflow
comparison pages, plus local validation and detection code.

## Proposed navigation

| Primary page | What a reader should accomplish | Content |
|---|---|---|
| Start here | Understand the tool and choose the next step | Short scope statement, supported input summary, links to run / inspect / assess, link to installation |
| Run ViralScan | Complete one recommended workflow | Preflight, reference prerequisites, combined host + virus build or reuse, quantification, validation, expected output paths |
| Inspect results | Find and read the main results | Report, per-virus summary, per-cell table, called-cell denominator, provenance receipts |
| Assess the evidence | Decide what a candidate supports | Technical integrity, molecule ambiguity, cell calling, read evidence, controls, limitations |
| Reference & help | Look up a specific detail | Installation alternatives, reference choices, commands, output columns, troubleshooting, advanced examples |

Use `index.md` for Start here and retain `quickstart.md` for Run ViralScan.
Add `inspect_results.md` and `assess_evidence.md`. Group existing detailed pages
under Reference & help with a shallow sidebar. Put developer API, changelog,
maintenance, migration, and support links on a secondary index or in the footer.
Keep notebooks reachable through one Examples page and links from relevant guides.

## Essential content

### Run ViralScan

1. Link to one recommended, verified installation route. State supported assay,
   chemistry, FASTQ pairing, and resource prerequisites before the command.
2. Show how to build a combined host + virus reference. Readers who already have
   a compatible index can skip that step. Link other reference modes separately.
3. Show one complete quantification command with consistent sample and directory
   names. Explain the called-cell requirement and how an external cell list fits.
4. Run `viralscan validate-run output/` and explain errors, warnings, and successful
   validation. State that it checks technical consistency, not biological validity.
5. Link directly to Inspect results with the expected report and table paths.

Keep cell-type enrichment, host response, alternate allocation methods, optional
alignment branches, and large molecule exports out of the first-run sequence.
Preserve the toy-data warning where relevant: `--cell-calling none` changes the
denominator and is not a routine workaround for failed cell calling.

### Inspect results

Start with a compact file / question / first-fields table. Cover `report.html`,
`results/viral_summary.tsv`, `results/per_cell_viral.tsv`, and
`results/called_cells.tsv`. Link AnnData and full schemas for downstream work.
Show one small, clearly labelled example row, including a weak or ambiguous case.
Use a verified fixture if available. Label invented teaching values as illustrative.

Explain molecule estimates versus reads, unique versus allocated support,
`detected` as a configured threshold result, and the called-cell denominator.
Explain `pct_infected_called` for within-run interpretation and the conditions
under which `pct_infected_comparable` is available for strategy comparisons.
Do not describe either field as proof of infection. Missing measurements mean
unavailable, not zero or a clean result.

### Assess the evidence

Use five checks, each with the relevant artifact, interpretation, and next action:

1. **Run integrity:** validation report, completion marker, input/reference/software
   provenance, warnings, and any provisional status. State what was actually checked.
2. **Cell support:** called-cell list, support in called cells versus background
   barcodes, chemistry/whitelist checks, and the configured detection thresholds.
3. **Mapping support:** unique and ambiguous molecule evidence, host competition,
   reference coverage, and relevant low-complexity or contamination diagnostics.
4. **Read inspection:** one verified `viralscan evidence` example, its exact
   reference/tool prerequisites, and the resulting competitive alignments and plots.
   Emphasize corrected molecules rather than PCR read multiplicity.
5. **Conclusion limits:** evidence tiers are rule-based categories, not calibrated
   probabilities; successful validation and threshold passage do not establish
   infection. Negative results depend on reference coverage and sensitivity.
   Link orthogonal controls and virus-specific caveats where they affect a decision.

Do not introduce new universal breadth, identity, or molecule cutoffs. Preserve
chemistry-specific and virus-specific limitations in the relevant reference pages.
Include sensitivity/control interpretation without promoting modelled bounds into
validated sensitivity claims.

## Content migration

| Existing material | Revision |
|---|---|
| Home page commands | Replace with the three reader tasks and concise scope |
| Installation | Show the recommended route first; group alternative environments below |
| Quickstart | Retain one complete main path; link alternative reference recipes |
| Output reference | Keep exact file/field definitions; move the reading sequence to Inspect results and evidence interpretation to Assess the evidence |
| FAQ | Retain symptom → cause → action troubleshooting; replace repeated explanations with links |
| Reference panel | Explain the general panel first; move extended anellovirus instructions to a clearly optional section |
| CLI reference | Preserve generated flag tables; move long study narratives out of command descriptions |
| Notebooks | Keep worked examples under one index; repair titles and distinguish synthetic illustrations from real-data validation |
| Contracts and maintenance | Retain as linked technical reference; remove internal task IDs and investigation history from the essential guides |

Keep one authoritative home for each explanation. Preserve scientific caveats and
provenance when shortening material. Retain existing page URLs and important anchors,
or supply redirects and update inbound links when splitting pages.

## Implementation sequence and acceptance

1. **Reconcile sources.** Identify the source and revision used for the published
   site. Its home page, theme, and commands differ from this local checkout, which
   also contains existing documentation changes. Preserve those edits and reconcile
   terminology against the intended release before rewriting examples.
2. **Restructure navigation.** Update the root toctree and home page; repair duplicate
   API entries and notebook titles. Keep the existing Sphinx stack unless it blocks
   the required navigation.
3. **Write the three core guides.** Complete Run, Inspect, and Assess with consistent
   examples and links. Readers must reach each directly from the home page.
4. **Trim and migrate.** Apply the content map, update links and anchors, and check
   destination pages for duplication. Do not hide necessary caveats in collapsed text.
5. **Verify the result.** Build Sphinx with warnings treated as errors; check links,
   generated CLI-reference consistency, and documentation contract tests. Validate
   examples against the target parser and an existing small fixture. A real FASTQ
   run needs suitable tools/data and is separate from a successful documentation
   build. Check mobile/sidebar readability and old deep links in the rendered site.

Done means a new reader can find the recommended run command, open the main
results, explain their denominator and units, distinguish a candidate from a
biological conclusion, and find the evidence needed for the next check without
reading API documentation or a notebook.

Rough effort: 2–3 working days for restructuring, editing, and validation, assuming
an existing working docs environment and no newly discovered command defects.
Publication is a separate step after reviewing the rendered changes.

## Review limits

Reviewed published HTML/text and local source/code. No visual browser or mobile
rendering check, workflow execution, notebook execution, or scientific performance
qualification was performed. This plan changes no commands, defaults, or results.
