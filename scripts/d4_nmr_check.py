# Naked mole-rat check: are the TCF12 / NBN / TCF4 retrocopy expansions real or assembly artifacts?
# Evidence: (1) counts in an independent NMR assembly (hetGla2) and in Damara mole-rat; (2) how the loci sit
# (scaffolds, spacing); (3) sequence: k-mer similarity between copies and between their 5 kb flanks
# (independent retro-insertions have unrelated flanks; segmental or assembly duplications share flanks), soft-masked share.
# Usage: d4_nmr_check.py <out_dir> <HLhetGla3.2bit> <genes comma> <control genes comma>
import csv, gzip, os, sys
from collections import defaultdict
from itertools import combinations
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from twobit import TwoBit

out_dir, twobit_path, genes_arg, control_arg = sys.argv[1:5]
symbols = genes_arg.split(",") + control_arg.split(",")
ASSEMBLIES = {"HLhetGla3": "Heterocephalus_glaber__naked_mole-rat__HLhetGla3",
	"hetGla2": "Heterocephalus_glaber__naked_mole-rat__hetGla2",
	"HLfukDam2": "Fukomys_damarensis__Damara_mole-rat__HLfukDam2",
	"cavPor3": "Cavia_porcellus__domestic_guinea_pig__cavPor3"}
FLANK = 5000
K = 15

gene_of_symbol = {}
for row in csv.DictReader(open(os.path.join(out_dir, "genes.tsv"), encoding="utf8"), delimiter="\t"):
	gene_of_symbol.setdefault(row["toga_symbol"], row["ensembl_gene"])

def read_loci(assembly_dir):
	path = os.path.join(out_dir, "pseudo_loci", assembly_dir + ".loci.tsv.gz")
	if not os.path.exists(path):
		# Assemblies outside the main species list (e.g. the second NMR assembly) live in pseudo_loci_extra
		path = os.path.join(out_dir, "pseudo_loci_extra", assembly_dir + ".loci.tsv.gz")
	if not os.path.exists(path):
		return None
	loci = defaultdict(list)
	with gzip.open(path, "rt", encoding="utf8") as f:
		next(f)
		for line in f:
			gene, chrom, start, end, strand, _ = line.rstrip("\n").split("\t")
			loci[gene].append((chrom, int(start), int(end), strand))
	return loci

def kmers(sequence):
	sequence = sequence.upper()
	return {sequence[k:k + K] for k in range(len(sequence) - K + 1) if "N" not in sequence[k:k + K]}

def jaccard(a, b):
	return len(a & b) / len(a | b) if a and b else 0.0

def revcomp(sequence):
	return sequence.upper().translate(str.maketrans("ACGTN", "TGCAN"))[::-1]

log = open(os.path.join(out_dir, "nmr_check.log"), "w", encoding="utf8")
def say(*items):
	text = " ".join(str(x) for x in items)
	print(text)
	log.write(text + "\n")

say("== Retrocopy loci per assembly ==")
loci_by_assembly = {label: read_loci(name) for label, name in ASSEMBLIES.items()}
for symbol in symbols:
	gene = gene_of_symbol[symbol]
	counts = []
	for label, loci in loci_by_assembly.items():
		counts.append(f"{label} {'NA' if loci is None else len(loci.get(gene, []))}")
	say(f"{symbol}: " + " | ".join(counts))

genome = TwoBit(twobit_path)
nmr = loci_by_assembly["HLhetGla3"]
say("\n== HLhetGla3: layout and sequence of the copies ==")
with open(os.path.join(out_dir, "nmr_check_loci.tsv"), "w", encoding="utf8") as table:
	table.write("symbol\tchrom\tstart\tend\tlength\tstrand\tscaffold_length\tsoftmasked_share\tgap_to_previous_same_gene\n")
	for symbol in symbols:
		gene = gene_of_symbol[symbol]
		items = sorted(nmr.get(gene, []))
		if not items:
			say(f"{symbol}: no loci")
			continue
		scaffolds = defaultdict(int)
		body_kmers = []
		flank_kmers = []
		gaps = []
		previous = None
		for chrom, start, end, strand in items:
			scaffolds[chrom] += 1
			body = genome.fetch(chrom, start, end)
			left = genome.fetch(chrom, start - FLANK, start)
			right = genome.fetch(chrom, end, end + FLANK)
			if strand == "-":
				body = revcomp(body)
				left, right = revcomp(right), revcomp(left)
			masked = sum(c.islower() for c in body) / max(len(body), 1)
			gap = start - previous[2] if previous and previous[0] == chrom else ""
			if gap != "":
				gaps.append(gap)
			previous = (chrom, start, end)
			body_kmers.append(kmers(body))
			flank_kmers.append(kmers(left) | kmers(right))
			table.write(f"{symbol}\t{chrom}\t{start}\t{end}\t{end - start}\t{strand}\t{genome.length(chrom)}\t{masked:.3f}\t{gap}\n")
		pairs = list(combinations(range(len(items)), 2))
		body_similarity = sorted(jaccard(body_kmers[a], body_kmers[b]) for a, b in pairs)
		flank_similarity = sorted(jaccard(flank_kmers[a], flank_kmers[b]) for a, b in pairs)
		lengths = sorted(end - start for _, start, end, _ in items)
		def summary(values):
			if not values:
				return "NA"
			return f"median {values[len(values) // 2]:.3f}, share >0.5: {sum(v > 0.5 for v in values) / len(values):.2f}, max {values[-1]:.3f}"
		say(f"{symbol}: {len(items)} loci on {len(scaffolds)} scaffolds (most on one scaffold: {max(scaffolds.values())}); "
			f"length median {lengths[len(lengths) // 2]} bp (min {lengths[0]}, max {lengths[-1]}); "
			f"same-scaffold gaps < 50 kb: {sum(g < 50000 for g in gaps)} of {len(gaps)}")
		say(f"   copy-vs-copy 15-mer Jaccard ({len(pairs)} pairs): {summary(body_similarity)}")
		say(f"   flank-vs-flank 15-mer Jaccard (5 kb each side): {summary(flank_similarity)}")
log.close()
