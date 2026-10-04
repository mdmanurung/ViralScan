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
    # HPV16 only (NC_001526); the "1618" in the prefix is a legacy file name.
    "HUM_PAP_1618": "Human papillomavirus 16",
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
    # Genomes added by the CAT-31/32/33 reconciliation (F-015, F-016). INSDC-only
    # records carry GenBank-style locus_tag gene ids (``AB027021_E1``), so the
    # accession itself is the panel token. Names use the clinically meaningful
    # papillomavirus type rather than the genus-level RefSeq ORGANISM string.
    "AB027021": "Human papillomavirus 82",
    "AB745694": "Human papillomavirus 160",
    # Token-prefixed gene ids from the same additions. These GTFs lead with a
    # virus/locus token rather than the accession (``OC43_EYW02_gp1``,
    # ``EBV2_HHV4tp2_gp1``, ``HIV1_gp1``), so the accession keys above do not match
    # them; the underscore-bounded prefix rule needs the token itself.
    "EBV2_HHV4tp2": "Human gammaherpesvirus 4 type 2",
    "FLUD_C1935": "Influenza D virus segment 7",
    "HBPV1_HBoV": "Bocaparvovirus primate 1",
    "HIV1": "Human immunodeficiency virus 1",
    "HIV2": "Human immunodeficiency virus 2",
    "HKU1": "Human coronavirus HKU1",
    "HMPV_D1Y22": "Human metapneumovirus",
    "HPV18": "human papillomavirus 18",
    "HPV31_QKA66": "human papillomavirus 31",
    "HPV33_QKB09": "human papillomavirus 33",
    "HPV35": "human papillomavirus 35",
    "HPV39": "human papillomavirus 39",
    "HPV45": "human papillomavirus 45",
    "HPV51": "human papillomavirus 51",
    "HPV52": "human papillomavirus 52",
    "HPV56": "human papillomavirus 56",
    "HPV58": "human papillomavirus 58",
    "HPV59": "human papillomavirus 59",
    "HPV66": "human papillomavirus 66",
    "HPV69": "human papillomavirus 69",
    "HPYV6_HPyV6": "Human polyomavirus 6",
    "HPYV7_HPyV7": "Human polyomavirus 7",
    "HTLV1": "Human T-cell leukemia virus type I",
    "HTLV2": "Human T-lymphotropic virus 2",
    "NL63": "Human coronavirus NL63",
    "OC43_EYW02": "Human coronavirus OC43",
    "OC43_YP": "Human coronavirus OC43",
    "SARS2_GU280": "SARS coronavirus 2",
    "SFV": "Western chimpanzee simian foamy virus",
    "TSPV_TSaPV": "Trichodysplasia spinulosa polyomavirus",
    # Accession-keyed ids from bundled HPV GTFs whose accessions had no map entry.
    # The underscore-bounded prefix rule then resolves ``NC_001355_HPV6bgp1`` etc.
    "NC_001355": "Human papillomavirus 6",
    "NC_001531": "Human papillomavirus 5",
    "NC_001576": "Human papillomavirus 10",
    "NC_001583": "Human papillomavirus 26",
    "NC_001586": "Human papillomavirus 32",
    "NC_001587": "Human papillomavirus 34",
    "NC_001591": "Human papillomavirus 49",
    "NC_001593": "Human papillomavirus 53",
    "NC_001595": "Human papillomavirus 7",
    "NC_001596": "Human papillomavirus 9",
    "NC_001676": "Human papillomavirus 54",
    "NC_001694": "Human papillomavirus 61",
    "NC_004104": "Human papillomavirus 90",
    "NC_004500": "Human papillomavirus 92",
    "NC_005134": "Human papillomavirus 96",
    "NC_038889": "Human papillomavirus 30",
    "NC_039089": "Human papillomavirus 71",
    "NC_075235": "Human papillomavirus 11",
    "NC_075251": "Human papillomavirus 38",
    "AF151983": "Human papillomavirus 83",
    "AF293960": "Human papillomavirus 84",
    "AF349909": "Human papillomavirus 86",
    "AF419318": "Human papillomavirus 91",
    "AF436128": "Human papillomavirus 89",
    "AF436130": "Human papillomavirus 74",
    "AJ400628": "Human papillomavirus 87",
    "AJ620205": "Human papillomavirus 43",
    "AJ620209": "Human papillomavirus 81",
    "AJ620211": "Human papillomavirus 94",
    "AY382778": "Human papillomavirus 93",
    "AY395706": "Human papillomavirus 62",
    "D21208": "Human papillomavirus 67",
    "DQ080079": "Human papillomavirus 68",
    "DQ080080": "Human papillomavirus 97",
    "DQ080082": "Human papillomavirus 106",
    "DQ080083": "Human papillomavirus 102",
    "EF422221": "Human papillomavirus 107",
    "EU410348": "Human papillomavirus 110",
    "EU410349": "Human papillomavirus 111",
    "FJ947080": "Human papillomavirus 115",
    "FM955837": "Human papillomavirus 98",
    "FM955838": "Human papillomavirus 99",
    "FM955839": "Human papillomavirus 100",
    "FM955840": "Human papillomavirus 104",
    "FM955841": "Human papillomavirus 105",
    "FM955842": "Human papillomavirus 113",
    "FN547152": "Human papillomavirus 125",
    "FN677755": "Human papillomavirus 150",
    "FN677756": "Human papillomavirus 151",
    "GQ244463": "Human papillomavirus 114",
    "GQ246950": "Human papillomavirus 117",
    "GQ246951": "Human papillomavirus 118",
    "GQ845442": "Human papillomavirus 120",
    "GQ845444": "Human papillomavirus 122",
    "GQ845446": "Human papillomavirus 124",
    "HE963025": "Human papillomavirus 159",
    "HF930491": "Human papillomavirus 174",
    "HM999995": "Human papillomavirus 143",
    "HM999997": "Human papillomavirus 145",
    "JF304768": "Human papillomavirus 152",
    "KC138720": "Human papillomavirus 78",
    "M12737": "Human papillomavirus 8",
    "M32305": "Human papillomavirus 47",
    "M73236": "Human papillomavirus 42",
    "U21941": "Human papillomavirus 70",
    "U31778": "Human papillomavirus 20",
    "U31779": "Human papillomavirus 21",
    "U31780": "Human papillomavirus 22",
    "U31781": "Human papillomavirus 23",
    "U31782": "Human papillomavirus 24",
    "U31783": "Human papillomavirus 28",
    "U31784": "Human papillomavirus 29",
    "U31785": "Human papillomavirus 36",
    "U31786": "Human papillomavirus 37",
    "U31788": "Human papillomavirus 44",
    "X55965": "Human papillomavirus 57",
    "X62843": "Human papillomavirus 13",
    "X74462": "Human papillomavirus 3",
    "X74466": "Human papillomavirus 12",
    "X74467": "Human papillomavirus 14",
    "X74468": "Human papillomavirus 15",
    "X74469": "Human papillomavirus 17",
    "X74470": "Human papillomavirus 19",
    "X74471": "Human papillomavirus 25",
    "X74473": "Human papillomavirus 27",
    "X74478": "Human papillomavirus 40",
    "X94164": "Human papillomavirus 72",
    "X94165": "Human papillomavirus 73",
    "Y15173": "Human papillomavirus 75",
    "Y15174": "Human papillomavirus 76",
    "Y15175": "Human papillomavirus 77",
    "Y15176": "Human papillomavirus 80",
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
    "HpV16": "Human papillomavirus 16",
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

# Viral families/genera with germline endogenous viral elements (EVEs) in the human
# genome, whose reads may come from host pre-mRNA rather than infection. Empty: no
# such family is in the shipped panel. Anelloviridae were listed until 2026-10-03
# (ANELLO-PRIOR.4) on an "EVEs in NALCN, LINC02742" basis that F-019 revised (>=90 %
# poly-G reads, <=10 % host-best); no germline human anellovirus EVE is known (one
# somatic integration in the SKNO-1 cell line, PMID 42671192). The mechanism stays
# for a real EVE family (e.g. inherited ciHHV-6, parked).
EVE_RISK_GENERA: frozenset[str] = frozenset()

# Families/genera whose calls carry a known read-artefact risk: low-complexity reads
# (poly-G no-signal reads, F-019; poly-A sinks on homopolymer-ended genomes, F-021)
# land on their low-complexity UTR/tail sequence. Reported as artifact_risk; a
# diagnostic label, never a filter (anelloviruses are commensal, ANELLO-PRIOR).
# ANELLO-PRIOR.3's measured read-level metric is meant to replace it.
LOW_COMPLEXITY_RISK_GENERA: frozenset[str] = frozenset(
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
