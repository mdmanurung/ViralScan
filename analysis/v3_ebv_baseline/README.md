# ViralScan 3.0 EBV molecule-counting baseline

This directory records the v3 `host-conservative` counting baseline for the
retained `SRR12682296` EBV run. The benchmark uses the production preparation
boundary (`bustools sort`, then streamed CB-UMI molecule resolution) and the
same allocator used by the CLI.

Final validation job: SLURM `25316336` (`COMPLETED`, exit `0:0`). The process
completed in 6 minutes 13 seconds. GNU time reported peak RSS of 5,644,620 KiB;
SLURM reported 5,278,448 KiB for the batch step. The allocator itself took
359.56 seconds.

The 103,145,071 BUS records contained 31,129,183 repeated read observations
that do not add molecule mass. The v3 audit reports:

- input molecules: 57,957,364
- unique molecules: 51,976,300
- ambiguous molecules: 4,220,847
- unresolved molecules: 1,760,217
- selected-matrix mass: 56,197,147

The invariants are exact:

`unique + ambiguous + unresolved = input`

`selected matrix = unique + allocated ambiguous = resolved`

For comparison only, the incompatible pre-v3 bustools matrix has mass
55,266,624. It is reported to demonstrate that the v3 selected matrix does not
reuse or numerically migrate legacy counts.

`host_conservative.json` is the machine-readable result. The resolved sorted
BUS and text files are retained for reproducibility. Job `25316225` is the
expected failed attempt that proved an unsorted raw BUS is rejected; job
`25316242` was an interim run before the final unique-layer independence fix.
