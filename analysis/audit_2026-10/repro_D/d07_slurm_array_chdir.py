"""write_slurm_array bakes a placeholder repo path into the job script."""
import os
import tempfile
from pathlib import Path

from viralscan.reference_strategy import write_slurm_array

p = write_slurm_array(Path(tempfile.mkdtemp()))
line = [l for l in p.read_text().splitlines() if "os.chdir" in l][0]
print(line, "| exists:", os.path.isdir(line.split('"')[1]))
