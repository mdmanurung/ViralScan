# Translated viral-protein preparation — 2026-10-14

User requested an attributed port of Laura Luebbert's kallisto translated
single-cell virus tutorial, then asked to include DNA viruses through their
expressed coding RNA. Scope: the local GTF corpus, all annotated viral CDS
proteins, separate from the production nucleotide detection workflow.

Implemented `scripts/prepare_viral_proteins.py`, the Snakemake workflow and
attribution under `analysis/translated_search/`, and four focused tests.
No aligner implementation, package dependency, production CLI change, sequence
download, commit, push, or real-sample analysis was added.

The local GTFs resolve to 393 catalogue-hash-matched cached GenBank records.
The package's minimal feature parser truncates multiline translations, so
the installed Biopython parser reads the full annotated translations instead.
Output: 4,124 proteins, 389 represented accessions, 296/299 catalogue labels;
492 proteins cover 111 RNA labels, 3,632 cover 185 DNA labels. Annotation-free
JN860200.1, D90457.1, M12737.1, and NC_005225.1 remain explicit coverage gaps;
modelled nucleotide features are not silently promoted to protein annotations.

Verification: four pytest passes; changed-file lint/format and data governance
pass. Fresh Snakemake preparation/index/count run succeeds with kallisto
0.52.0, kb-python 0.30.2, and Snakemake 9.23.1. Synthetic reads encoding Aichi
and HSV-2 proteins recover exactly 20 UMIs each, with no other positive virus
counts, using a synthetic host mask. Receipts and failed-attempt logs remain
under the ignored `analysis/translated_search/work/smoke_2026-10-14/` root.
kb-python 0.30.2 lacks `kb count -k`; the index carries k=31. Its informational
CLI commands return status 1, so version provenance uses installed metadata
from the configured Python environment. Snakemake-created empty output parents
are accepted; nonempty reference directories are preserved.

Prepared a populated ignored config using available human-only genome/cDNA
inputs. Its production DAG passes dry-run; full human-mask indexing and
real-sample runs are not executed. Protein matches concern expressed coding
RNA; noncoding latency RNAs and reactivation claims require other evidence.
