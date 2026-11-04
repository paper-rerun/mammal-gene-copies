# Mammal gene copy number, retrocopies and gene loss vs. human (TOGA), joined to lifespan and body mass

Version 1.0.0, built on 2026-11-04. Part of the Paper Rerun project (https://paper-rerun.github.io). Licence: CC BY 4.0 (see `LICENSE.md`); please also cite the sources listed under "How to cite".

This dataset puts three things TOGA already computes, but stores in separate files, into one table: for every human (hg38) protein-coding gene that TOGA projects (19,467 genes) and each of 562 mammal assemblies (one per species), it gives the TOGA gene status, how many orthologous loci TOGA found, and how many processed-pseudogene (retrocopy) loci TOGA annotated. It adds AnAge life-history traits, the matching tip of a public mammal tree, two assembly-quality flags, and a per-assembly outlier ranking.

It was built to look at Peto's paradox (large, long-lived mammals without more cancer). Read "What this table can't show" before using it for that.

## Files

| File | Rows | What it is |
|---|---|---|
| `out/gene_by_species.tsv.gz` | 10,940,454 (562 x 19,467) | long table, one row per assembly x gene |
| `out/gene_by_species.parquet` | same | same table, typed, `NA` stored as null |
| `out/species.tsv` | 562 | one row per assembly: taxonomy, traits, tree tip, QC |
| `out/genes.tsv` | 19,467 | one row per reference gene: symbols and gene-set membership |
| `out/outliers_by_assembly.tsv` | top 50 per assembly and metric | per-lineage outlier ranking |
| `out/pseudo_loci/*.loci.tsv.gz` | 488 files | merged retrocopy loci with coordinates. **Not in this repository** (0.66 GB); `scripts/d1_pseudogene_loci.py` rebuilds them, and `out/pseudo_loci_md5.txt` lets you check the result |
| `out/nmr_check*.{log,tsv}` | | naked mole-rat artifact check (see below) |
| `figures/` | | five summary figures, see `figures/README.md` |
| `scripts/` | | everything needed to rebuild |

### `gene_by_species.tsv.gz`

| Column | Meaning |
|---|---|
| `assembly` | TOGA assembly name; join key to `species.tsv` |
| `ensembl_gene` | human reference gene (Ensembl ID as used by TOGA); join key to `genes.tsv` |
| `toga_status` | TOGA gene classification from `loss_summ_data.tsv.gz`: `I` intact, `PI` partially intact, `UL` uncertain loss, `L` lost, `M` missing, `PM` partially missing, `PG` paralogous projection only, `NA` not classified |
| `orthology_class` | set of TOGA orthology classes over the gene's transcripts, `|`-separated: `one2one`, `one2many`, `many2one`, `many2many`, `one2zero` (no intact ortholog) |
| `n_ortholog_loci` | number of distinct query loci TOGA reports as orthologous (distinct `q_gene` in `orthologsClassification.tsv.gz`, all transcripts). 0 for `one2zero`. Values of 2+ are candidate duplications, but see the false-duplication flag |
| `n_retrocopy_loci` | number of processed-pseudogene loci TOGA annotated for this gene: projections of all its transcripts from `processedPseudogeneAnnotation.bed.gz`, overlapping intervals merged. `NA` when TOGA has no such file for the assembly (74 assemblies, HTTP 404) |

### `species.tsv`

`assembly`, `species`, `common_name`, `toga_order_folder` (TOGA download folder: an order or a clade such as Whippomorpha), `ncbi_accession_or_source`, `contig_n50`, `scaffold_n50` (both from TOGA `overview.table.tsv`), `anage_name`, `max_longevity_yrs`, `adult_weight_g` (AnAge adult weight, or AnAge body mass when adult weight is empty), `anage_data_quality`, `anage_sample_size`, `upham2019_tip` (tip label in the Upham et al. 2019 tree, empty if not matched), and:

| Column | Meaning |
|---|---|
| `has_retrocopy_scan` | 1 if TOGA's processed-pseudogene file exists for this assembly |
| `n_retrocopy_loci_total` | sum of `n_retrocopy_loci` over all genes. A locus annotated against several genes of one family is counted once per gene, so this is not a genome-wide retrocopy count |
| `qc_retrocopy_scan_sparse` | 1 if `n_retrocopy_loci_total` < 10,000. The log10 histogram of this total has a separate cluster at 10^2.75 to 10^3.75 (98 assemblies) against a main mode at 10^5 to 10^5.5. It is not clade biology (it splits genera) and not file date. Retrocopy counts from flagged assemblies are not comparable with the rest |
| `frac_multicopy` | share of genes with at least one ortholog that have `n_ortholog_loci` >= 2 |
| `frac_lost` | share of `L` among genes classified `I`, `PI` or `L` |
| `frac_missing` | share of genes classified `M` or `PM` |
| `qc_false_duplication` | 1 if `frac_multicopy` > 3 x the median over all 562 assemblies (median 7.11%, threshold 21.3%): 25 assemblies. Such assemblies most likely carry uncollapsed haplotypes; their 2+ ortholog loci are mostly not real duplications. `frac_multicopy` does not track contig N50 (Spearman 0.008 in the pilot), so N50 alone does not catch them |
| `qc_false_duplication_pilot35` | the pilot's looser list: 35 assemblies with more than 14.5% of an 890-gene panel multi-copy (3 x that panel's median of 4.8%). All 25 strict-flag assemblies are in it |

### `genes.tsv`

`ensembl_gene`, `toga_symbol` (symbol embedded in the TOGA transcript name), `hgnc_symbol`, `hgnc_locus_group` (current HGNC), and 0/1 gene-set columns:

| Column | Source | Genes in the table |
|---|---|---|
| `set_dna_repair_wood` | Wood, Lange and Wood human DNA repair genes (MD Anderson), updated Nov 2025 | 218 |
| `set_tp53_reactome` | Reactome R-HSA-3700989 "Transcriptional Regulation by TP53", Reactome release 98 (CC0) | 358 |
| `set_tumour_suppressor_intogen` | IntOGen Compendium Cancer Genes, release 2024-09-20, genes whose only role is LoF (CC0) | 166 |
| `set_genage_human` | GenAge human, build 21 | 304 |

The pilot used KEGG hsa04115 for the TP53 pathway. It was replaced here because of KEGG's redistribution terms.

### `outliers_by_assembly.tsv`

For each assembly and each metric (`n_retrocopy_loci`, `n_ortholog_loci`), the 50 genes where this assembly sits furthest above the others. `robust_z = (value - median_others) / (1.4826 x max(MAD_others, 0.5))`. "Others" are all other assemblies, except that assemblies flagged `qc_retrocopy_scan_sparse` are left out of the retrocopy baseline and assemblies flagged `qc_false_duplication` are left out of the ortholog-loci baseline. Flagged assemblies are still ranked, and the flags are repeated in each row. `n_others_at_or_above` and `n_same_order_at_or_above` tell a one-assembly event from a clade-wide one (both elephants, all cetaceans).

Positive control: TP53 in the African elephant has 22 retrocopy loci against a baseline median of 1, and no non-proboscidean assembly has more than 10. Genome-wide it ranks 13th of 19,467 genes in the African elephant and 35th in the Asian elephant (`out/check_table.log`). Genes shared by both elephants rank above it (TNFAIP8L3 with 33 and 35 loci, MSRB1 with 21 and 24). These are candidate proboscidean retrocopy expansions that were not checked further.

### Naked mole-rat check (`out/nmr_check.log`, `out/nmr_check_loci.tsv`, `out/nmr_check_hetGla2_layout.log`)

The table's naked mole-rat assembly is hetGla2 (highest contig N50). There TCF12 (49 loci), NBN (33) and TCF4 (27) are among its top retrocopy outliers. Checked against the second, independent assembly (HLhetGla3, DNA Zoo, chromosome-length scaffolds), the Damara mole-rat and the guinea pig:

| Gene | hetGla2 | HLhetGla3 | Damara mole-rat | guinea pig | hetGla2 locus length, median | Verdict |
|---|---|---|---|---|---|---|
| NBN | 33 | 9 | 0 | 1 | 2,471 bp | looks real |
| TCF12 | 49 | 2 | 2 | 0 | 340 bp | likely artifact or repeat-derived fragments |
| TCF4 | 27 | 0 | 1 | 1 | 679 bp | likely artifact or repeat-derived fragments |
| GAPDH (control) | 32 | 11 | 23 | 34 | 1,111 bp | known retrocopied gene |

HLhetGla3 has a sparse retrocopy scan (4,409 loci in total against 240,577 in hetGla2), so its counts are lower for every gene, including the controls. NBN still reaches 9 there, close to GAPDH (11). In HLhetGla3 the 9 NBN copies are about full length (1,654 to 2,515 bp, median 2,442), spread over 8 chromosome-length scaffolds of 50 to 140 Mb, with 0 of 1 same-scaffold pairs closer than 50 kb. They are not repeat-masked (0 to 3.1% soft-masked), and their 5 kb flanks are unrelated (15-mer Jaccard median 0.000, max 0.007). That is the signature of independent retro-insertions, not of a segmental or assembly duplication. The copies differ from each other (copy-vs-copy Jaccard median 0.166), so they are not recent. TCF12 and TCF4 show up only in hetGla2, mostly as short fragments, which fits an annotation or repeat artifact better than an expansion. Nothing here says whether the NBN copies are expressed or functional.

## Sources and versions

| Source | Version | Licence | md5 |
|---|---|---|---|
| TOGA, human hg38 reference (Kirilenko et al. 2023, Science 380:eabn3107), genome.senckenberg.de/download/TOGA/human_hg38_reference | per-assembly files dated 2021-11-19 to 2025-04-18; `overview.table.tsv` | no formal licence; the TOGA authors confirmed reuse with citation (Kirilenko et al. 2023, plus the download link) in https://github.com/hillerlab/TOGA/issues/226 (2026-11-04) | `overview.table.tsv` 55617f38853abc427c779f3f8f700d72; the 1,124 `orthologsClassification` and `loss_summ_data` files: `out/toga_tables_md5.txt`. The processed-pseudogene files were streamed and not stored; md5 of the derived loci files: `out/pseudo_loci_md5.txt` |
| AnAge (HAGR) | build 15, 2023-07-03 | HAGR conditions of use | `anage_dataset.zip` 7ec4d19f6c8f57b22595a9ffd50ea66b |
| GenAge human | build 21, 2023-08-28 | HAGR conditions of use | `genage_human.zip` 2172ffe8a8f85557c8f9d27cf4eb2be5 |
| IntOGen drivers | 2024-09-20 | CC0 | `IntOGen-Drivers-20240920.zip` e0d7aaf8f9195dc2e3e6d9784f325570 |
| HGNC complete set | downloaded 2026-11-04 | CC0 | `hgnc_complete_set.txt` 39d46ee3ced2310056ec206019b67425 |
| Wood DNA repair list | page updated Nov 2025, downloaded 2026-11-04 | public web page | `wood_dna_repair.html` 7017bbc3b5f80ab664b6f3be5d76fea3 |
| Reactome pathways (GMT) | release 98 | data CC0 | `ReactomePathways.gmt.zip` 1f695ad09eb3ad5d0331ddddf7c2372d |
| Upham, Esselstyn and Jetz 2019 mammal tree, DNA-only 4,098 species, node-dated, MCC | `MamPhy_fullPosterior_BDvr_DNAonly_4098sp_topoFree_NDexp_MCC_v2_target.tre` | see VertLife | bdbc71bf37afa3dee72d536f424a6da1 |
| Naked mole-rat assembly HLhetGla3 (TOGA DNA Zoo 2bit, only for the check) | 2022-11-15 | DNA Zoo | 1c3d61e3b1fa96eb4dcbc634ca39e80e |

All md5 values are also in `out/source_md5.tsv`. Assembly choice: one assembly per species, the one with the highest contig N50 in `overview.table.tsv`. 21 newer primate T2T folders missing from that table were skipped.

## How to rebuild

Python 3.12. Scripts without third-party packages run as `py -3.12 -I -X utf8`. `d3`, `d5` and `d6` need numpy, pyarrow and matplotlib (`.venv`). Put the downloads in `dl/`: `scripts/s2_select_assemblies.py` picks one assembly per species from TOGA's `overview.table.tsv`, and `scripts/s3_fetch_toga.sh out/toga_species.tsv dl` fetches the TOGA tables. AnAge, the gene sets and the tree come from the URLs under "Sources and versions" (md5 in `out/source_md5.tsv`). These download steps were run in an earlier pilot folder and were not re-run from this layout. Then, from this folder:

```
py -3.12 -I -X utf8 scripts/d0_transcript_map.py out/toga_species.tsv dl/toga out/reference_transcripts.tsv
py -3.12 -I -X utf8 scripts/d1_pseudogene_loci.py out/toga_species.tsv out/reference_transcripts.tsv out/pseudo_loci dl/toga
py -3.12 -I -X utf8 scripts/d2_build.py dl dl out
py -3.12 -I -X utf8 scripts/d2b_pilot_flag.py out/species.tsv out/exclude_high_multicopy.txt
.venv/Scripts/python.exe -I scripts/d3_outliers.py out 50
.venv/Scripts/python.exe -I scripts/d5_parquet.py out
.venv/Scripts/python.exe -I scripts/check_table.py out
.venv/Scripts/python.exe -I scripts/d6_figures.py out dl/hagr/anage_data.txt figures
```

`d1` streams about 20 GB of TOGA pseudogene files with 6 processes and keeps only merged loci (about 0.6 GB). `check_table.py` must print "all checks passed" after any rebuild. It checks row counts, the elephant TP53 retrocopies, UCP1, F12 and KLKB1 lost in all 9 cetaceans checked (F12 and KLKB1 intact in hippo, cow and dog), UCP1 lost in pig, MMP20 lost in three baleen whales and intact in two toothed whales, and TP53 in the outlier ranking (top 20 in the African elephant, top 50 in the Asian elephant). The naked mole-rat check is `d4_nmr_check.py` and `d4b_nmr_layout_hetgla2.py`.

## How to cite

Please cite the sources this table is built from:

- TOGA: Kirilenko B. M. et al. (2023). Integrating gene annotation with orthology inference at scale. Science 380, eabn3107. Data: https://genome.senckenberg.de/download/TOGA/human_hg38_reference
- AnAge and GenAge (HAGR, CC BY 3.0): de Magalhães J. P. et al. (2024). Human Ageing Genomic Resources: updates on key databases in ageing research. Nucleic Acids Research 52(D1), D900-D908.
- Mammal tree: Upham N. S., Esselstyn J. A., Jetz W. (2019). Inferring the mammal tree. PLOS Biology 17(12), e3000494.
- IntOGen Compendium of Cancer Genes (CC0), Reactome release 98 (CC0), HGNC (CC0).

## What this table can't show

- **Correlation across species misses lineage events.** In the pilot, a phylogenetic regression of TP53 retrocopy count on body mass over 303 species gave t = -0.28. The elephant TP53 expansion happened once, in one lineage, so it does not move a regression over hundreds of species. In the same pilot, DNA-repair, TP53-pathway and tumour-suppressor genes tracked mass or lifespan no better than a random control set, and the single "significant" set-level result was a great-ape clade effect. Use the outlier ranking, or clade-aware tests, for this kind of question.
- **False duplications.** Some assemblies carry uncollapsed haplotypes, and up to 41.6% of their genes look multi-copy (`frac_multicopy`). Use `qc_false_duplication` before reading anything into `n_ortholog_loci` >= 2.
- **Retrocopy scans are not uniform.** 98 assemblies have very sparse processed-pseudogene annotations (`qc_retrocopy_scan_sparse`), and 74 have none. The two naked mole-rat assemblies differ by a factor of about 55 in total retrocopy loci (240,577 in hetGla2 against 4,409 in HLhetGla3). Compare retrocopy counts only between non-sparse assemblies, and preferably within a clade.
- **No sequence, no function.** TOGA gives coordinates for retrocopies, not sequence or reading frame. A retrocopy here can be a full-length gene, a dead pseudogene or a fragment. Protein abundance and expression, the basis of e.g. the bowhead CIRBP result (Firsanov et al. 2025), do not appear in any genome annotation.
- **Human-referenced.** Genes absent from human are invisible. TOGA cannot tell a loss in the query from a gain in the human lineage without an outgroup (see the TOGA README).
- **TOGA quality limits.** Gene status depends on assembly contiguity: fragmented assemblies produce more `M` and spurious `L` calls (Spearman between contig N50 and lost fraction -0.215 in the pilot). `UL` should be read as unknown, not as lost.
- **Traits.** AnAge covers about two thirds of the assemblies, with uneven data quality. Maximum longevity depends on how many animals were observed.
- **Everything here is correlational.** A copy-number difference that lines up with lifespan is a lead for experiments, not a cause.
