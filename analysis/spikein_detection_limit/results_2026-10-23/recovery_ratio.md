# Flank / old net kb recovery ratio

Ratio = flank mean / old mean of per-seed net kb recovery. RATIO_MIN = 0.8; flagged only where the Q15 verdict is keep and the old recovery is >= 0.1. The keep/retire verdict is never changed here.

| background | layout | level | old mean | flank mean | ratio | seed min | seed max | Q15 | needs user decision | NA reason |
|---|---|---|---|---|---|---|---|---|---|---|
| D03 | c005 | 3 | 0.556 | 0.556 | 1 | 1 | 1 | retire |  |  |
| D03 | c005 | 10 | 0.567 | 0.783 | 1.38 | 1.17 | 1.58 | retire |  |  |
| D03 | c005 | 30 | 0.685 | 0.772 | 1.13 | 1.1 | 1.15 | retire |  |  |
| D03 | c005 | 100 | 0.649 | 0.791 | 1.22 | 1.15 | 1.3 | retire |  |  |
| D03 | c005 | 300 | 0.641 | 0.787 | 1.23 | 1.21 | 1.25 | retire |  |  |
| D03 | c100 | 3 | 0.556 | 0.556 | 1 | 1 | 1 | retire |  |  |
| D03 | c100 | 10 | 0.567 | 0.783 | 1.38 | 1.17 | 1.58 | retire |  |  |
| D03 | c100 | 30 | 0.685 | 0.772 | 1.13 | 1.1 | 1.15 | retire |  |  |
| D03 | c100 | 100 | 0.649 | 0.778 | 1.2 | 1.16 | 1.24 | retire |  |  |
| D03 | c100 | 300 | 0.641 | 0.784 | 1.22 | 1.2 | 1.25 | retire |  |  |
| R01 | c005 | 3 | 0.889 | 1 | 1.12 | 1 | 1.5 | retire |  |  |
| R01 | c005 | 10 | 0.867 | 1 | 1.15 | 1.11 | 1.25 | retire |  |  |
| R01 | c005 | 30 | 0.778 | 0.956 | 1.23 | 1.17 | 1.32 | retire |  |  |
| R01 | c005 | 100 | 0.83 | 0.97 | 1.17 | 1.13 | 1.22 | retire |  |  |
| R01 | c005 | 300 | 0.843 | 0.961 | 1.14 | 1.1 | 1.17 | retire |  |  |
| R01 | c100 | 3 | 0.889 | 1 | 1.12 | 1 | 1.5 | retire |  |  |
| R01 | c100 | 10 | 0.867 | 1 | 1.15 | 1.11 | 1.25 | retire |  |  |
| R01 | c100 | 30 | 0.778 | 0.956 | 1.23 | 1.17 | 1.32 | retire |  |  |
| R01 | c100 | 100 | 0.833 | 0.973 | 1.17 | 1.13 | 1.21 | retire |  |  |
| R01 | c100 | 300 | 0.843 | 0.961 | 1.14 | 1.1 | 1.17 | retire |  |  |

kb molecules per spiked read, net of the unspiked replicate, uncapped; seeds are pseudo-replicates of one source library
