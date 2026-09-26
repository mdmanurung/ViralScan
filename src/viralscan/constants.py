"""Shared constants for ViralScan.

The ``VIRUS_NAME_MAP`` dictionary maps the abbreviated prefix used in the
bundled GTF gene IDs to the human-readable virus name. It used to be
duplicated verbatim in ``detection.py`` and ``umap.py``; both modules now
import it from here.
"""

VIRUS_NAME_MAP: dict[str, str] = {
    "AICHI": "Aichi virus",
    "AUSBATLYSSA": "Australian Bat Lyssavirus",
    "BANNA": "Banna virus",
    "BARMAH": "Barmah forest virus",
    "BKPOLY": "BK polyomavirus",
    "BUNYAMW": "Bunyamwera virus",
    "BUNYA": "Bunyavirus La Crosse",
    "CERC_HERP": "Cercopithecine herpesvirus",
    "CHIKUNG": "Chikungunya virus",
    "COSA_A": "Cosavirus A",
    "COWPOX": "Cowpox virus",
    "COXSACKIE": "Coxsackievirus",
    "CRIMEAN": "Crimean Congo hemorrhagic fever virus",
    "EQUINE_ENCE": "Eastern_equine_encephalitis_virus",
    "EBOLA": "Ebolavirus",
    "ECHO": "Echovirus",
    "ENCEPHAL": "Encephalomyocarditis virus",
    "EPSTEIN": "Epstein-Barr virus",
    "EURBATLYSSA": "European bat lyssavirus",
    "GB": "GB virus C_Hepatitis G virus",
    "HANTAAN": "Hantaan virus",
    "HENDRA": "Hendra virus",
    "HEP_A": "Hepatitis A virus",
    "HEP_B": "Hepatitis B virus",
    "HEP_C": "Hepatitis C virus",
    "HEP_DELTA": "Hepatitis delta virus",
    "HEP_E": "Hepatitis E virus",
    "HUM_ADENO": "Human adenovirus",
    "HUM_ASTRO": "Human astrovirus",
    "HUM_COR": "Human coronavirus",
    "HUM_CYTO": "Human cytomegalovirus",
    "HUM_ENTERO": "Human enterovirus 68, 70",
    "HUM_HERP1": "Human herpesvirus 1",
    "HUM_HERP2": "Human herpesvirus 2",
    "HUM_HERP6B": "Human herpesvirus 6b",
    "HUM_HERP6": "Human herpesvirus 6",
    "HUM_HERP7": "Human herpesvirus 7",
    "HUM_HERP8": "Human herpesvirus 8",
    "HUM_PAP_1618": "Human papillomavirus 16,18",
    "HUM_PAP_1": "Human papillomavirus 1",
    "HUM_PAP_2": "Human papillomavirus 2",
    "HUM_PARA": "Human parainfluenza",
    "HUM_PARVO": "Human parvovirus B19",
    "HUM_RESP": "Human respiratory syncytial virus",
    "HUM_RHINO": "Human rhinovirus",
    "HUM_SARS": "Human SARS coronavirus",
    "INFL_A": "Influenza A virus",
    "INFL_B": "Influenza B virus",
    "INFL_C": "Influenza C virus",
    "JAP_ENCE": "Japanese encephalitis virus",
    "JC_POLY": "JC polyomavirus",
    "KI_POLY": "KI Polyomavirus",
    "LAKE_VIC": "Lake Victoria marburgvirus",
    "LANGAT": "Langat virus",
    "LASSA": "Lassa virus",
    "LOUPING": "Louping ill virus",
    "LYMPH": "Lymphocytic choriomeningitis virus",
    "MAYARO": "Mayaro virus",
    "MEASLES": "Measles virus",
    "MERKEL": "Merkel cell polyomavirus",
    "MERS": "MERS coronavirus",
    "MOLLU": "Molluscum contagiosum virus",
    "MONKEYPOX": "Monkeypox virus",
    "MUMPS": "Mumps virus",
    "MUR_VAL": "Murray valley encephalitis virus",
    "NIPAH": "Nipah virus",
    "NORWALK": "Norwalk virus",
    "ORF": "Orf virus",
    "OROPOU": "Oropouche virus",
    "ONYONG": "O'nyong-nyong virus",
    "POLIO": "Poliovirus",
    "RABIES": "Rabies virus",
    "ROSA_A": "Rosavirus A",
    "ROSS_RIVER": "Ross river virus",
    "ROTA_A": "Rotavirus A",
    "ROTA_B": "Rotavirus B",
    "ROTA_C": "Rotavirus C",
    "RUBELLA": "Rubella virus",
    "SALI_A": "Salivirus A",
    "SAPPORO": "Sapporo virus",
    "SEMLIKI": "Semliki forest virus",
    "SEOUL": "Seoul virus",
    "SINDBIS": "Sindbis virus",
    "ST_LOUIS": "St. louis encephalitis virus",
    "TICK": "Tick-borne powassan virus",
    "TTV": "Torque teno virus",
    # Anelloviridae genera — added when the expanded anellovirus reference was
    # incorporated (clareaulab/anellovirus_reference, 2025).  These entries
    # cover accession-keyed gene IDs (e.g. "NC_014076.2_gene1") that the
    # boundary-aware prefix rule in virus_grouping resolves via the full
    # accession key in anellovirus.anello_name_map() — the entries here serve
    # as human-readable labels for any prefix-based fallback if needed, and are
    # also exposed as the canonical genus display names.
    "Alphatorquevirus": "Alphatorquevirus",
    "Betatorquevirus": "Betatorquevirus",
    "Gammatorquevirus": "Gammatorquevirus",
    "Samektorquevirus": "Samektorquevirus",
    "Memtorquevirus": "Memtorquevirus",
    "Hetorquevirus": "Hetorquevirus",
    "Gyrovirus": "Gyrovirus",
    "Anelloviridae": "Anelloviridae",
    "TOSCANA": "Toscana virus",
    "VACCINIA": "Vaccinia virus",
    "VARICELLA": "Varicella-zoster virus",
    "VARIOLA": "Variola virus",
    "VEN_EQU": "Venezuelan equine encephalitis virus",
    "VES_STOM": "Vesicular stomatitis virus",
    "WES_EQU": "Western equine encephalitis virus",
    "WES_NILE": "West Nile virus",
    "WU_POLY": "WU polyomavirus",
    "YABA": "Yaba-like disease virus",
    "YELLOW": "Yellow fever virus",
    "ZIKA": "Zika virus",
    # ── Alternate panel naming schemes ──────────────────────────────────────
    # The Serratus-derived panels shipped alongside the bundled GTFs use a
    # different token convention from the bundled VIRUS_NAME_MAP keys. Measured
    # against the real `log/analysis.txt` of the merged EBV/HSV-1/HHV-6B
    # benchmark and the covid PBMC runs, 215/4845 (4.4%) of viral gene IDs
    # resolved to *no* virus name, and the covid run published 9 of 17 rows of
    # `viral_summary.tsv` under raw gene IDs (`HHV1gp00p39`, `CeHV2gUL24`,
    # `MPXV_gp132`). Each unresolved gene became its own "virus", so
    # per-virus aggregation, `accession_breadth`, sibling cross-mapping and
    # `eve_risk` were all computed per gene instead of per virus.
    #
    # These keys are the underscore-bounded leading tokens of those schemes.
    # Every one of them is separated from the gene token by `_`, so they match
    # under the strict boundary rule in `virus_grouping` without weakening it.
    "ADENO": "Adeno-associated virus",  # panel AAV2; distinct from HUM_ADENO
    "DENV": "Dengue virus",
    "DUGBE": "Dugbe virus",
    "IMMUNO": "Human immunodeficiency virus",  # panel token: IMMUNO_HIV1*
    "JUNIN": "Junin virus",
    "MACH": "Machupo virus",
    "MOKO": "Moko virus",
    "PICHI": "Pichivirus",
    "PUUM": "Puumala virus",
    "RIFT": "Rift Valley fever virus",
    "SIMI": "Simian foamy virus",
    "TLYMPHO": "Human T-lymphotropic virus",  # panel token: TLYMPHO_HTLV1*
    "UUKU": "TTV-like mini virus",
    "YABAM": "Yaba-like disease virus",  # panel token: YABAM_YMTV*
    # Human herpesviruses under the ICTV HHV-n convention, used by the
    # STARsolo-packaged viral panel. Only the underscore-bounded forms land
    # here; the concatenated forms (``HHV5wtgp045``) are aliases below.
    "HHV3": "Varicella-zoster virus",
    "HHV7": "Human herpesvirus 7",
    "HHV8": "Human herpesvirus 8",
    "HAdVC": "Human adenovirus",
    "MARV": "Lake Victoria marburgvirus",
    "MCPyV": "Merkel cell polyomavirus",
    "B19V": "Human parvovirus B19",
    "EBLV1": "European bat lyssavirus",
    "WUPyV": "WU polyomavirus",
    "NoVGI": "Norwalk virus",
    "RVA": "Rotavirus A",
}

#: Gene-ID prefixes for panel schemes that concatenate the virus token and the
#: gene token with **no separator** (``Ydvgp129``, ``TTVgp1``, ``HHV1gp00p39``)
#: or that use a token with an embedded strain/segment suffix
#: (``HHV5wtgp045``, ``HHV8GK18_gp56``, ``VACWR202``).
#:
#: These cannot live in :data:`VIRUS_NAME_MAP`: its boundary rule deliberately
#: rejects a key followed by a letter, and that rule is load-bearing. It is what
#: stops ``EPSTEIN_HHV4_BORF1`` being mislabelled Orf virus (key ``ORF`` hidden
#: inside ``BORF1``) and ``BUNYAMW_...`` being read as Bunyavirus La Crosse.
#: Relaxing it re-opens both bugs.
#:
#: Instead these are matched by plain prefix, longest key first, and **only
#: after** the strict boundary rule has matched nothing. The escape hatch is
#: therefore strictly additive: no gene ID that resolves today can change name.
#:
#: Only tokens with an unambiguous standard assignment are listed. Panel tokens
#: whose virus cannot be identified with confidence are deliberately absent and
#: keep the raw-gene-ID fallback, which at least stays visible in the output:
#: ``QKL08``, ``HRCV``, ``KPV``, ``G128``, ``SCV12``, ``unassigned``.
#:
#: Sorted longest-first at match time, so a longer sibling (``HHV6B``) is
#: preferred over a shorter one (``HHV6``) without needing a special case.
VIRUS_GENE_ID_ALIASES: dict[str, str] = {
    "AAV2": "Adeno-associated virus",
    "CeHV2": "Cercopithecine herpesvirus",
    "FLUAV": "Influenza A virus",
    "FLUBV": "Influenza B virus",
    "FLUCV": "Influenza C virus",
    "HCoV229E": "Human coronavirus 229E",
    "HHV1": "Human herpesvirus 1",
    "HHV2": "Human herpesvirus 2",
    "HHV4": "Epstein-Barr virus",
    "HHV5": "Human cytomegalovirus",
    "HHV6": "Human herpesvirus 6",
    "HHV6B": "Human herpesvirus 6b",
    "HHV7": "Human herpesvirus 7",
    "HHV8": "Human herpesvirus 8",
    "Hpv1": "Human papillomavirus 1",
    "HpV16": "Human papillomavirus 16,18",
    "HpV2": "Human papillomavirus 2",
    "IMMUNO": "Human immunodeficiency virus",
    "MOCV": "Molluscum contagiosum virus",
    "MPXV": "Monkeypox virus",
    # Non-anellovirus accession-keyed gene IDs (`{acc}_geneN`, emitted by
    # ncbi_fetch._genome_to_gtf). Anellovirus accessions are resolved in tier 1
    # by anellovirus.anello_name_map(); this is where other accessions belong.
    "NC_045512": "SARS coronavirus 2",
    "TLYMPHO": "Human T-lymphotropic virus",
    "TTV": "Torque teno virus",
    "VACW": "Vaccinia virus",
    "VARV": "Varicella-zoster virus",
    "YdV": "Yaba-like disease virus",
    "Ydv": "Yaba-like disease virus",  # panel token is lower-case v
}

# Ensembl species registry used by `viralscan build-ref`.
# Keys are the short names accepted on the CLI (case-insensitive).
# Values are (ensembl_species_name, genome_assembly_name) tuples.
# Species names match the directory layout at
#   https://ftp.ensembl.org/pub/current_fasta/<species_name>/cdna/
ENSEMBL_SPECIES: dict[str, tuple[str, str]] = {
    "human": ("homo_sapiens", "GRCh38"),
    "mouse": ("mus_musculus", "GRCm39"),
    "rat": ("rattus_norvegicus", "mRatBN7.2"),
    "zebrafish": ("danio_rerio", "GRCz11"),
    "chicken": ("gallus_gallus", "bGalGal1"),
    "macaque": ("macaca_mulatta", "Mmul_10"),
    "pig": ("sus_scrofa", "Sscrofa11.1"),
    "dog": ("canis_lupus_familiaris", "ROS_Cfam_1.0"),
    "ferret": ("mustela_putorius_furo", "MusPutFur1.0"),
    "cat": ("felis_catus", "Felis_catus_9.0"),
    "marmoset": ("callithrix_jacchus", "mCalJac1.pat.X"),
    "cow": ("bos_taurus", "ARS-UCD1.3"),
    "sheep": ("ovis_aries", "ARS-UI_Ramb_v2.0"),
    "rabbit": ("oryctolagus_cuniculus", "OryCun2.0"),
    "hamster": ("mesocricetus_auratus", "MesAur1.0"),
    "drosophila": ("drosophila_melanogaster", "BDGP6.46"),
    "celegans": ("caenorhabditis_elegans", "WBcel235"),
    "xenopus": ("xenopus_tropicalis", "UCB_Xtro_10.0"),
}

# Sibling virus pairs that share sufficient nucleotide identity to cause
# reciprocal multimapper allocation under the global EM. HHV-6A/6B share
# ~95% identity; HSV-1/2 share ~80%. When one sibling dominates by
# SIBLING_CROSSMAP_RATIO_THRESHOLD or more, the weaker signal is likely
# EM bleed from shared-region multimappers rather than genuine co-infection.
#
# NOT the same list as `sibling_virus_pairs` in
# analysis/v3_validation/protocol.yaml, despite the shared name. That one is an
# evaluation population for endpoints D13/D24/E5 and deliberately spans a
# relatedness gradient, including EBV/KSHV — two gammaherpesviruses in different
# genera, chosen precisely because they are the most distant pair. Adding EBV/KSHV
# here would be wrong: they do not cross-map, so the ratio heuristic would
# annotate genuine EBV+KSHV co-infection (common in KS and PEL) as EM bleed.
# Membership here requires near-identity; membership there requires only that the
# pair be prespecified.
SIBLING_VIRUS_PAIRS: dict[str, str] = {
    "Human herpesvirus 6": "Human herpesvirus 6b",
    "Human herpesvirus 6b": "Human herpesvirus 6",
    "Human herpesvirus 1": "Human herpesvirus 2",
    "Human herpesvirus 2": "Human herpesvirus 1",
}

SIBLING_CROSSMAP_RATIO_THRESHOLD: float = 50.0

# Viral families/genera known to have endogenous viral element (EVE) integrations
# in the human genome. Reads mapping to these viruses may originate from intronic
# pre-mRNA of expressed host genes (e.g. Anellovirus EVEs in NALCN, LINC02742)
# rather than exogenous infection, especially when counts concentrate on 1-2 loci.
EVE_RISK_GENERA: frozenset[str] = frozenset(
    {
        "Alphatorquevirus",
        "Betatorquevirus",
        "Gammatorquevirus",
        "Samektorquevirus",
        "Memtorquevirus",
        "Hetorquevirus",
        "Gyrovirus",
        "Anelloviridae",
        "Torque teno virus",
    }
)
