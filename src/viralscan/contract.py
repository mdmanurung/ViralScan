"""Version constants of ViralScan's public contract, in one place.

Each is part of what a run writes or what a reader accepts. A change here is a contract
change: bump it deliberately and add a migration note (docs/migration.md, docs/release_policy.md).
"""

# schema_version of every JSON document ViralScan writes (run_manifest.json, run_complete.json,
# the validate-run report, the evidence manifest, the positive-control receipt).
SCHEMA_VERSION = "3.0.0"
# count_schema_version stamped into adata.uns of the multimap h5ad; readers refuse any other value.
COUNT_SCHEMA_VERSION = "3.0.0"
# Output layout: the version written, and the versions readers still accept.
LAYOUT_VERSION = 3
SUPPORTED_LAYOUT_VERSIONS = (1, 2, 3)
# results/artifact_inventory.json.
INVENTORY_SCHEMA_VERSION = "1.1.0"
