# ViralScan v3 truth-panel status

**Status: not yet run.** The authoritative tasks are `SCI-01` through `SCI-05`
and `VAL-01` through `VAL-10` in [`../../PLAN.md`](../../PLAN.md). Freeze the
schema-valid protocol and its digest before examining new v3 outcomes.

The active machine-readable draft and its explicit training/holdout gates are in
[`../v3_validation/protocol.yaml`](../v3_validation/protocol.yaml), with usage in
[`../v3_validation/README.md`](../v3_validation/README.md).

Prior releases demonstrated workflows on published public scRNA-seq datasets,
but those results are not v3 validation and do **not** estimate formal
false-positive or false-negative rates against known single-cell infection
truth. Until the v3 truth panel is complete, the README and manuscript must
state that FPR/FNR are not formally characterized.

The minimum legacy scaffold below is retained for context; it is not the v3
acceptance contract:

1. A negative-control host-only scRNA-seq dataset processed with the same host+virus
   reference (measures the single-UMI false-positive rate).
2. A synthetic spike-in or read-mixing dataset with known viral read/barcode
   assignments (measures recovery / false-negative rate at controlled abundance).
3. A positive-control infected dataset with an orthogonal published infection label
   or an accepted threshold.

The v3 panel additionally requires preregistered train/holdout partitions,
supported chemistry geometries, sibling-virus and host-homology adversarial
conditions, deterministic molecule truth, sample-level uncertainty, and the
complete harmonized workflow matrix defined in `PLAN.md`.
