from viralscan.scripts.ncbi_fetch import (
    NONCODING_FEATURE_KEY,
    NCBIFetchError,
    _gene_token,
    _gtf_attributes,
    _locus_fields,
    _noncoding_token,
    _panel_gene_ids,
    _parse_location,
    _transcript_order,
    iter_features,
)


def _genbank_to_gtf(genbank_text: str, accession: str) -> str:
    """Minimal GenBank → GTF converter.

    Extracts CDS features and, because the most abundant viral transcripts are
    routinely non-coding, ``misc_RNA`` features, and emits one ``exon`` line per
    location interval with ``gene_id`` and ``transcript_id`` attributes. This is
    the minimum that ``kb ref``/``gffread`` need to extract transcript sequences.

    ``misc_RNA`` is the key EBV annotates EBER1/EBER2 under — RefSeq
    ``NC_007605.1`` 6629..6795 and 6956..7128, ``NC_009334.1`` 6634..6800 and
    6961..7133 — and those are the highest-abundance latent EBV transcripts, so
    a CDS-only reference cannot see latent infection at all.  ``ncRNA``/miRNA
    features are deliberately *not* emitted: EBV's are 18–24 nt, shorter than
    kallisto's ``k=31``, so they contain no 31-mer and would be dead targets.

    CDS rows are emitted first in flatfile order and the non-coding rows after
    them, so every CDS gene ID, transcript ID, row order and byte position is
    identical to a CDS-only build.  That is what keeps the 1,994-accession
    Anelloviridae catalogue and every already-cached ``<accession>.gtf``
    untouched; it also matches how the vendor EBV reference is assembled, with
    the two EBER records appended after the 94 CDS records.

    A single-exon non-coding feature is emitted as the one exon NCBI annotates;
    no second exon is fabricated.  Single-exon targets are not at risk from
    ``kb ref``: the vendor EBV reference already carries both EBERs as
    single-exon ``exon`` records, 2,444 of the 2,511 packaged Anelloviridae
    genes are single-exon, and the single-exon whole-genome placeholder
    ``MW455439.1_gene1`` took 1,167,103 UMI in a COVID run.

    Two consequences of admitting the whole ``misc_RNA`` class are worth
    knowing, both measured across the 2,249 cached flatfiles (6 records carry
    any ``misc_RNA`` at all; the 1,994-accession Anelloviridae panel carries
    none, so its GTFs and the packaged catalogue are unaffected):

    * ``k`` is the length floor, not the exon count.  HHV-8 ``NC_009333.1``
      annotates 10 miRNAs as ``misc_RNA`` at 20–23 nt, which contain no 31-mer,
      so ``kallisto index`` skips them; they stay in the GTF so the annotation
      records that NCBI calls them, rather than hiding a filter in the writer.
    * Two records annotate a genome-scale ``misc_RNA``: hepatitis A
      ``NC_001489.1`` 1..7478 is the whole 7,478-nt genome, and HTLV-2
      ``NC_001488.1`` 316..8751 is 94 % of its 8,952 nt.  Emitted into a
      ``viral_gene`` index these are exactly the whole-genome pseudo-transcripts
      that collapsed 99.8 % of anellovirus UMI into one bucket, so resolving
      them belongs to the two-index-tier split in
      ``docs/plans/2026-09-27-viral-reference-panel-expansion.md`` §5.1, not to
      a span threshold invented here.

    Gene IDs are genome-scoped — see the module docstring.  Raises
    :class:`NCBIFetchError` when the record carries neither a CDS nor a
    ``misc_RNA`` feature, which is the caller's signal to fall back to
    :func:`_whole_genome_gtf_from_fasta`.
    """
    locus = _locus_fields(genbank_text)
    seqid_field = str(locus["version"]) or accession
    genome_length = int(locus["genome_length"])
    cds_features = [
        (location, qualifiers)
        for key, location, qualifiers in iter_features(genbank_text)
        if key == "CDS"
    ]
    noncoding_features = [
        (location, qualifiers)
        for key, location, qualifiers in iter_features(genbank_text)
        if key == NONCODING_FEATURE_KEY
    ]
    if not cds_features and not noncoding_features:
        raise NCBIFetchError(
            f"No CDS or {NONCODING_FEATURE_KEY} features found in GenBank record for "
            f"{accession}; cannot build a GTF for kb ref."
        )

    out: list[str] = []
    seen_tokens: dict[str, int] = {}
    for group, token_for, biotype in (
        (cds_features, _gene_token, "protein_coding"),
        (noncoding_features, _noncoding_token, NONCODING_FEATURE_KEY),
    ):
        for gene_id, transcript_id, token, location, qualifiers in _panel_gene_ids(
            group, accession, token_for=token_for, seen=seen_tokens
        ):
            exons = _transcript_order(_parse_location(location))
            if not exons:
                continue
            attrs = _gtf_attributes(
                gene_id, transcript_id, token, qualifiers, exons, genome_length, biotype
            )
            for start, end, strand in exons:
                out.append(f"{seqid_field}\tNCBI\texon\t{start}\t{end}\t.\t{strand}\t0\t{attrs}")

    if not out:
        raise NCBIFetchError(
            f"CDS and {NONCODING_FEATURE_KEY} features in GenBank record for {accession} "
            "carry no parseable location; cannot build a GTF for kb ref."
        )
    return "\n".join(out) + "\n"
