# Reference transcript -> reference gene map (TOGA hg38 reference), union over all assemblies' orthologsClassification
# Usage: d0_transcript_map.py <toga_species.tsv> <toga_download_dir> <out.tsv>
import csv, gzip, os, sys

species_path, toga_dir, out_path = sys.argv[1:4]
gene_by_transcript = {}
for row in csv.DictReader(open(species_path, encoding="utf8"), delimiter="\t"):
	path = os.path.join(toga_dir, os.path.basename(row["toga_dir"]), "orthologsClassification.tsv.gz")
	with gzip.open(path, "rt", encoding="utf8") as f:
		next(f)
		for line in f:
			gene, transcript = line.split("\t", 2)[:2]
			gene_by_transcript[transcript] = gene

with open(out_path, "w", encoding="utf8") as out:
	out.write("t_transcript\tt_gene\tsymbol\n")
	for transcript in sorted(gene_by_transcript):
		out.write(f"{transcript}\t{gene_by_transcript[transcript]}\t{transcript.split('.', 1)[1]}\n")
print("transcripts", len(gene_by_transcript), "genes", len(set(gene_by_transcript.values())))
