# Truth Panel Plan for ViralScan Specificity/Sensitivity (not yet run)

The current release validates against published public scRNA-seq datasets but does
**not** estimate formal false-positive or false-negative rates against a known
single-cell infection truth. Until a truth panel is complete, the README and the
manuscript must state that FPR/FNR are not formally characterized.

Minimum acceptable truth panel:

1. A negative-control host-only scRNA-seq dataset processed with the same host+virus
   reference (measures the single-UMI false-positive rate).
2. A synthetic spike-in or read-mixing dataset with known viral read/barcode
   assignments (measures recovery / false-negative rate at controlled abundance).
3. A positive-control infected dataset with an orthogonal published infection label
   or an accepted threshold.

Status: planning only — not part of the current validation package.
