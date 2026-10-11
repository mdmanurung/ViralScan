"""Structural checks on constants.VIRUS_NAME_MAP / VIRUS_GENE_ID_ALIASES."""
import ast
import collections
from pathlib import Path

import viralscan.constants as c

src = Path(c.__file__).read_text()
tree = ast.parse(src)
for node in tree.body:
    target = getattr(node, "target", None) or (node.targets[0] if isinstance(node, ast.Assign) else None)
    if isinstance(node, (ast.AnnAssign, ast.Assign)) and isinstance(node.value, ast.Dict):
        keys = [k.value for k in node.value.keys if isinstance(k, ast.Constant)]
        dup = [k for k, n in collections.Counter(keys).items() if n > 1]
        print(target.id, "n_keys", len(keys), "duplicate literal keys:", dup)

# One virus under several display names (same casefolded, whitespace-normalised value)
by_norm = collections.defaultdict(set)
for k, v in c.VIRUS_NAME_MAP.items():
    by_norm[v.casefold().replace("_", " ")].add(v)
print("same name, different spelling:", {k: sorted(v) for k, v in by_norm.items() if len(v) > 1})

# Same biological virus, different map values
groups = {
    "HIV-1": ["IMMUNO", "HIV1"],
    "HTLV-1": ["TLYMPHO", "HTLV1"],
    "HPV 16": ["HUM_PAP_1618", "HpV16"],
}
for g, keys in groups.items():
    vals = {}
    for k in keys:
        vals[k] = c.VIRUS_NAME_MAP.get(k) or c.VIRUS_GENE_ID_ALIASES.get(k)
    print(g, vals, "-> one name?" , len(set(vals.values())) == 1)

# keys of the two maps with different values
print("conflicts NAME_MAP vs ALIASES:",
      {k: (c.VIRUS_NAME_MAP[k], v) for k, v in c.VIRUS_GENE_ID_ALIASES.items()
       if k in c.VIRUS_NAME_MAP and c.VIRUS_NAME_MAP[k] != v})
