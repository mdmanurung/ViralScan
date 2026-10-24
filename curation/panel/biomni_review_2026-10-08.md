# External review of the panel-expansion plan (Biomni, 2026-10-08) and what we do about it

Source: Biomni task `tsk_015noA6X3BFfq7SsTft6FMKj` (project "ViralScan — viral reference completeness...", model tier
`max`; the underlying model cannot be pinned). It is external model output with literature citations that were not
re-checked here; the dispositions below are ours. "Verified" means checked against this repo's code, not against the
literature.

| # | Sev. | Finding (our words) | Disposition |
|---|---|---|---|
| 1 | high | The plan never says how multi-gene UMIs become calls; that decides whether cross-mapping loses or inflates counts. | Partly covered: the package pins `DEFAULT_MULTIMAP_METHOD = "host-conservative"` (`defaults.py:5`). **Gap confirmed:** the planned gates measure pseudoalignment, not the production counting path. Adopt: run the separability check through the real counting command and report unique-UMI and EC-level evidence. |
| 2 | high | A host D-list in a combined index may mask host targets. | Not on the default path (docs prefer `--host-filter starsolo`; the D-list measured weakest). Applies only to `--genome-dlist` builds: adopt a host-gene count-pinning gate if that path is ever the production one. |
| 3 | high | No redundancy cap outside anelloviruses; close relatives turn unique UMIs into ambiguous ones and lower counts for existing viruses. | Adopt. Per-family identity cap with named exceptions (dengue serotypes, HCV genotypes, HPV types, HIV-1 groups, HBV genotypes; influenza one per subtype/lineage; one SARS-CoV-2 reference, per `CAT-21`/`CAT-22`). |
| 4 | high | Host cross-talk is measured but not budgeted; no count of host-virus mixed equivalence classes. | Adopt as a WP3 gate (target about 0 with a documented exception list). |
| 5 | med | Error-free fragments overstate separability. | Adopt: simulate reads with errors, 3' bias and poly-A, per segment. |
| 6 | med | md5 twin detection misses rotated or reverse-complemented circular genomes. | Adopt for circular families (canonical rotation + strand before hashing). Current `panel_candidates.py` checks exact sequences only; `twin_checked` says which rows it could check. |
| 7 | med | Circular genomes lose reads spanning the origin in a linear index. | **Not found** in the build code (no wrap bases appended; only feature annotation handles origin-spanning joins). Check, then append the first k-1 bases for circular records if confirmed. |
| 8 | med | No genome quality gate for GenBank-only candidates. | Adopt: prefer RefSeq; near-complete vs the species exemplar; host/vector screen. |
| 9 | low | The EC-health gate needs definition; `--ec-max-size` defaults to no maximum. | Adopt: also track total EC count, index size, RAM, tool versions, index checksum; define what the existing "discarded ECs" figure counts. |
| 10 | low | Non-polyadenylated genomes (flaviviruses) are under-captured by oligo-dT. | Adopt as documentation: expected capture per family. |
| 11 | high | Host-field relevance (H1) is noisy; do not auto-promote. | Already how `panel_candidates.py` works (no row is promoted automatically). Extend: record the evidence source per accepted row. |
| 12 | high | The pool has holes (herpesviruses, polyomaviruses, HBV/HDV, HIV/HTLV, enteroviruses, HAV/HEV, parvoviruses, adenoviruses, poxviruses, paramyxo-/pneumoviruses, arboviruses, pegiviruses, circovirus, redondoviruses, gut viruses). | Adopt: build the pool top-down and diff against it. The census already shows gaps (human parvovirus 4, polyomavirus 9, Astrovirus VA1, Zika `NC_035889`, HCV genotype 7). |
| 13 | high | "Exclude HERV/EVE risk" is not operational, and a blanket rule would drop bornaviruses. | Adopt: pre-agree an action rule (mask / demote to decoy / exclude) on the host-overlap measurement; never add HERVs as viruses; document ciHHV-6 as an interpretation caveat. `risk_class=eve` is used by nothing yet. |
| 14 | med | Use ICTV VMR as the spine, joined to Virus-Host DB and the NCBI census; Serratus only to nominate. | Adopt (needs a download; ICTV VMR and Virus-Host DB are not in the repo). |
| 15 | low | Deprioritise giant viruses and comprehensive phages. | Adopt. |
| 16 | med | The vector/reagent set is missing phiX174 (and MS2, Ad5, baculovirus, vaccinia/MVA, VSV, Sabin polio, YF-17D, WPRE/SV40 cassettes). | Adopt: define an explicit reagent set, one row per element with accession.version. |
| 17 | med | `role` must be machine-actionable (infection-candidate, reagent, vector, vaccine-strain, commensal, decoy, uncertain); lentiviral vector vs HIV-1 must report as ambiguous. | Adopt. **Verified:** `role` is carried into `results/virus_identity.tsv` (`virus_identity.py:575`) but nothing else reads it, and only the 9 murine retroviruses use `decoy`. |
| 18 | high | WP2 curates before separability is measured. | **Adopt, changes the order:** run a k-mer / pilot-index separability check on the candidate superset between WP1 and WP2. |
| 19 | high | The freeze gate depends on thresholds and homology actions that are deferred. | Adopt: freeze criteria need interim thresholds and a stated homology action rule, or calibration moves ahead of the freeze. |
| 20 | high | Silent-change watchlist: masking existing genomes, twin consolidation, added relatives lowering counts, NCBI version drift, sibling-group edits, tool versions, D-list. | Adopt: frozen manifest (accession.version, md5, role, sibling group, tool versions, index checksum) plus a signed before/after count-diff on a pinned dataset. |
| 21 | med | A list of acceptance tests for "safe to freeze". | Adopt as the WP5 checklist (reconciliation, manifest, realistic self-consistency, mixed-EC budget, host-only negatives real and simulated, 2-3 public positives, EC health, count diff, role docs, named sign-off). |
| 22 | low | Pin the census date and inputs; decide the census question first. | Done: `census_raw.json` carries the query and date; the census question was answered "yes" before this review. |

Revised order: WP1 (done) -> **WP1b separability on the candidate superset** -> WP2 curation (top-down pool, roles,
caps) -> WP3 gates -> WP4 tests/docs/packaging -> WP5 user-owned build and freeze.
