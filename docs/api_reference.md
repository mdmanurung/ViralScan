# API reference (autodoc)

Generated from the NumPy-style docstrings of the importable modules. The
hand-written [API guide](api.md) has worked examples; this page is the
signature-level reference. Snakemake worker scripts (`detection`, `multimap`,
`umap`, ...) are workflow internals, not a stable API.

## Cell-type enrichment

```{eval-rst}
.. autofunction:: viralscan.enrichment.cell_type_enrichment
.. autofunction:: viralscan.enrichment.write_cell_type_enrichment
```

## Reference data cache

```{eval-rst}
.. autofunction:: viralscan.data_fetch.fetch_viral_data
.. autofunction:: viralscan.data_fetch.ensure_viral_data
.. autofunction:: viralscan.data_fetch.cache_valid
```

## Reference building

```{eval-rst}
.. autofunction:: viralscan.scripts.ncbi_fetch.fetch_reference
.. autofunction:: viralscan.scripts.build_reference.fetch_host_cdna
```
