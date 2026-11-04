"""Draft source-licence records for the reference manifest (PLAN REF-02 / REF-05).

Every entry is a DRAFT summary of the source's published terms. None is legal advice and none
has been reviewed by a maintainer: ``status`` stays ``"unreviewed"`` and ``review_required``
stays true until a maintainer signs off. The table in ``docs/reference_licences.md`` carries
the same sources with more detail; keep the two in step.
"""

from typing import Any

DRAFT_PREFIX = "DRAFT, needs maintainer/legal review: "

#: source label (the manifest ``source`` field) -> (terms summary, URL of the published terms)
_DRAFTS: dict[str, tuple[str, str]] = {
    "NCBI nuccore": (
        "NCBI places no restrictions on use of GenBank/RefSeq sequence data; a submitter may "
        "claim patent or copyright on a record, so redistribution is not cleared per record.",
        "https://www.ncbi.nlm.nih.gov/home/about/policies/",
    ),
    "Ensembl cDNA": (
        "Ensembl data are freely available without restriction under the EMBL-EBI terms of use; "
        "cite the Ensembl release.",
        "https://www.ensembl.org/info/about/legal/disclaimer.html",
    ),
    "local_input": (
        "User-supplied FASTA/GTF; its terms are unknown to ViralScan and are the user's to check.",
        "",
    ),
}


def draft_source_licence(source: str) -> dict[str, Any]:
    """Return the manifest ``source_licence`` record for *source*, still unreviewed."""
    terms, url = _DRAFTS.get(source, ("", ""))
    text = f"{DRAFT_PREFIX}{terms} {url}".strip() if terms else None
    return {"status": "unreviewed", "terms": text, "source": source, "review_required": True}
