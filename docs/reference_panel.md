# Reference Panel

ViralScan uses **195 viral GTF annotation files** distributed through Zenodo
under DOI `10.5281/zenodo.20112332`. Fetch the panel once after installation:

```bash
viralscan data fetch
```

The command verifies the Zenodo checksum and unpacks the GTF files to
`~/.cache/viralscan/data/`. These annotations were generated from RefSeq
GenBank entries using the `Viral_GTF_maker.py` script
(whole-genome-as-gene model: one gene per chromosome, covering the full
replicon length).

---

## Included viruses

The panel covers the following virus families and genera (non-exhaustive):

- **Adenoviridae** — Human adenovirus
- **Anelloviridae** — Torque teno virus (Alphatorquevirus homin spp.)
- **Arenaviridae** — Lassa virus, Lymphocytic choriomeningitis virus
- **Bunyaviridae / Hantaviridae** — Hantaan virus, Seoul virus, Bunyamwera virus, La Crosse virus
- **Caliciviridae** — Norwalk virus, Sapporo virus
- **Coronaviridae** — SARS-CoV, MERS-CoV, Human coronavirus 229E / HKU1 / NL63 / OC43
- **Filoviridae** — Ebolavirus, Lake Victoria marburgvirus
- **Flaviviridae** — Dengue virus, Hepatitis C virus, West Nile virus, Yellow fever virus, Zika virus
- **Herpesviridae** — EBV, CMV, HSV-1/2, HHV-6/7/8, VZV
- **Orthomyxoviridae** — Influenza A, B, C
- **Papillomaviridae** — Human papillomavirus 1, 2, 16, 18
- **Paramyxoviridae** — Measles virus, Mumps virus, Nipah virus, Hendra virus
- **Parvoviridae** — Human parvovirus B19
- **Picornaviridae** — Poliovirus, Rhinovirus, Coxsackievirus, Echovirus, Enterovirus 68/70
- **Polyomaviridae** — BK polyomavirus, JC polyomavirus, Merkel cell polyomavirus
- **Poxviridae** — Vaccinia virus, Cowpox virus, Monkeypox virus, Variola virus
- **Reoviridae** — Rotavirus A/B/C, Banna virus
- **Rhabdoviridae** — Rabies virus, Australian bat lyssavirus
- **Togaviridae** — Chikungunya virus, Rubella virus, Alphavirus spp.
- **Other** — Hepatitis A/B/D/E, Adeno-associated virus (AAV), and more

---

## Adding custom references

You can supplement the cached panel with your own annotation files using the
`-gtf` flag (comma-separated list) together with `--reference`:

```bash
viralscan --reference \
  -fasta custom_virus.fasta \
  -gtf   custom_virus.gtf \
  -o output/ \
  -s1 R1.fastq.gz -s2 R2.fastq.gz
```

Alternatively, use `--ncbi-accession` to fetch and build a reference for any
RefSeq nucleotide accession on-the-fly.

---

## GTF format

Each panel GTF file follows the RefSeq GTF convention:

```
NC_001477.1  RefSeq  gene  1  10735  .  +  .  gene_id "DENV_DV1_gp1"; ...
```

The `gene_id` prefix (e.g. `DENV`) is used by `detection.py` to look up the
human-readable virus name in `VIRUS_NAME_MAP` (`constants.py`).

---

## Cache location

By default, `viralscan data fetch` writes to `~/.cache/viralscan/data/`.
Use `viralscan data fetch --cache-dir PATH` to populate another cache root.
