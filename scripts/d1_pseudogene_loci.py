# Stream every species' TOGA processedPseudogeneAnnotation.bed.gz and keep merged retrocopy loci per reference gene.
# Projections of all transcripts of one gene that overlap on the same strand-agnostic interval are one locus.
# Output per assembly: <out_dir>/<assembly_dir>.loci.tsv.gz with gene, symbol, chrom, start, end, strand, n_projections
# Usage: d1_pseudogene_loci.py <toga_species.tsv> <transcript_map.tsv> <out_dir> [local_toga_dir]
import csv, gzip, os, sys, urllib.request
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor

BASE = "https://genome.senckenberg.de/download/TOGA/human_hg38_reference/"
species_path, map_path, out_dir = sys.argv[1:4]
local_dir = sys.argv[4] if len(sys.argv) > 4 else None
os.makedirs(out_dir, exist_ok=True)

gene_by_transcript = {}
for row in csv.DictReader(open(map_path, encoding="utf8"), delimiter="\t"):
	gene_by_transcript[row["t_transcript"]] = row["t_gene"]

def open_stream(toga_dir):
	name = os.path.basename(toga_dir)
	if local_dir:
		local_path = os.path.join(local_dir, name, "processedPseudogeneAnnotation.bed.gz")
		if os.path.exists(local_path):
			return gzip.open(local_path, "rt", encoding="utf8")
	request = urllib.request.urlopen(BASE + toga_dir + "/processedPseudogeneAnnotation.bed.gz", timeout=300)
	return gzip.open(request, "rt", encoding="utf8")

def process(row):
	name = os.path.basename(row["toga_dir"])
	out_path = os.path.join(out_dir, name + ".loci.tsv.gz")
	if os.path.exists(out_path):
		return name, "cached"
	intervals = defaultdict(list)
	try:
		with open_stream(row["toga_dir"]) as stream:
			for line in stream:
				fields = line.split("\t")
				transcript = fields[3].rsplit(".", 1)[0]
				gene = gene_by_transcript.get(transcript, "NA:" + transcript.split(".", 1)[1])
				intervals[gene].append((fields[0], int(fields[1]), int(fields[2]), fields[5].strip()))
	except Exception as error:
		return name, "FAIL " + repr(error)
	n_loci = 0
	with gzip.open(out_path + ".part", "wt", encoding="utf8") as out:
		out.write("ensembl_gene\tchrom\tstart\tend\tstrand\tn_projections\n")
		for gene in sorted(intervals):
			items = sorted(intervals[gene])
			merged = []
			for chrom, start, end, strand in items:
				if merged and merged[-1][0] == chrom and start <= merged[-1][2]:
					merged[-1][2] = max(merged[-1][2], end)
					merged[-1][4] += 1
				else:
					merged.append([chrom, start, end, strand, 1])
			for chrom, start, end, strand, count in merged:
				out.write(f"{gene}\t{chrom}\t{start}\t{end}\t{strand}\t{count}\n")
			n_loci += len(merged)
	os.replace(out_path + ".part", out_path)
	return name, n_loci

rows = list(csv.DictReader(open(species_path, encoding="utf8"), delimiter="\t"))
if __name__ == "__main__":
	with ProcessPoolExecutor(max_workers=6) as pool:
		for name, result in pool.map(process, rows):
			print(name, result, flush=True)
