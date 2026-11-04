# Validate a (re)built dataset against known biology. Exits with an error on the first failed check.
# Usage: check_table.py <out_dir>
import csv, gzip, os, sys
from collections import defaultdict

out_dir = sys.argv[1]
TP53, UCP1, F12, KLKB1, MMP20 = "ENSG00000141510", "ENSG00000109424", "ENSG00000131187", "ENSG00000164344", "ENSG00000137674"
CETACEANS = ["HLbalMys1", "HLeubGla1", "balAcu1", "HLbalMus2", "HLmegNov2", "HLphyCat2", "HLorcOrc2", "HLturTru4", "HLdelLeu2"]
PROBOSCIDEA = {"HLloxAfr4", "HLeleMax1"}
WANTED = {TP53, UCP1, F12, KLKB1, MMP20}

species = list(csv.DictReader(open(os.path.join(out_dir, "species.tsv"), encoding="utf8"), delimiter="\t"))
genes = list(csv.DictReader(open(os.path.join(out_dir, "genes.tsv"), encoding="utf8"), delimiter="\t"))
assemblies = {s["assembly"] for s in species}
dog = [s["assembly"] for s in species if s["species"] == "Canis lupus familiaris"]
NON_CETACEAN_INTACT = ["HLhipAmp3", "bosTau9"] + dog

cell = {}
rows_per_assembly = defaultdict(int)
n_rows = 0
with gzip.open(os.path.join(out_dir, "gene_by_species.tsv.gz"), "rt", encoding="utf8") as f:
	header = next(f).rstrip("\n").split("\t")
	assert header == ["assembly", "ensembl_gene", "toga_status", "orthology_class", "n_ortholog_loci", "n_retrocopy_loci"], header
	for line in f:
		assembly, gene, status, orthology, n_loci, n_retro = line.rstrip("\n").split("\t")
		n_rows += 1
		rows_per_assembly[assembly] += 1
		if gene in WANTED:
			cell[(assembly, gene)] = (status, n_loci, n_retro)

def check(condition, message):
	if not condition:
		sys.exit("FAILED: " + message)
	print("ok:", message)

check(len(species) == len(assemblies), f"species table has {len(species)} unique assemblies")
check(n_rows == len(species) * len(genes), f"long table rows {n_rows} = {len(species)} assemblies x {len(genes)} genes")
check(all(v == len(genes) for v in rows_per_assembly.values()), "every assembly has one row per gene")

# Positive control: elephant TP53 retrocopies (Sulak et al. 2016 report 19 retrogenes in loxAfr3)
lox = int(cell[("HLloxAfr4", TP53)][2])
ele = int(cell[("HLeleMax1", TP53)][2])
others = [int(v[2]) for (a, g), v in cell.items() if g == TP53 and a not in PROBOSCIDEA and v[2] != "NA"]
check(lox >= 19, f"African elephant TP53 retrocopy loci = {lox} (>= 19 published; 22 in the 2026-11-04 build)")
check(ele >= 14, f"Asian elephant TP53 retrocopy loci = {ele} (>= 14 published)")
check(max(others) < min(lox, ele), f"no non-proboscidean assembly reaches the elephants (max {max(others)})")
check(cell[("HLloxAfr4", TP53)][0] in ("I", "PI"), "elephant TP53 parent gene intact")

# Known losses
for gene, name in ((UCP1, "UCP1"), (F12, "F12"), (KLKB1, "KLKB1")):
	statuses = {a: cell[(a, gene)][0] for a in CETACEANS}
	check(all(s == "L" for s in statuses.values()), f"{name} lost (L) in all 9 cetaceans {statuses}")
for gene, name in ((F12, "F12"), (KLKB1, "KLKB1")):
	statuses = {a: cell[(a, gene)][0] for a in NON_CETACEAN_INTACT}
	check(all(s in ("I", "PI") for s in statuses.values()), f"{name} intact in hippo, cow, dog {statuses}")
check(cell[("susScr11", UCP1)][0] == "L", "UCP1 lost in pig (known loss)")
mysticetes = {a: cell[(a, MMP20)][0] for a in ("balAcu1", "HLbalMus2", "HLeubGla1")}
check(all(s == "L" for s in mysticetes.values()), f"MMP20 lost in minke, blue and right whale {mysticetes}")
odontocetes = {a: cell[(a, MMP20)][0] for a in ("HLorcOrc2", "HLturTru4")}
check(all(s in ("I", "PI") for s in odontocetes.values()), f"MMP20 intact in orca and dolphin {odontocetes}")

# Outlier ranking finds the positive control. Genome-wide (19,467 genes) the 2026-11-04 build ranks TP53 13th in the
# African elephant and 35th in the Asian elephant; other proboscidean retrocopy families (TNFAIP8L3, PCGF5, MSRB1) rank above it
outliers_path = os.path.join(out_dir, "outliers_by_assembly.tsv")
if os.path.exists(outliers_path):
	ranks = {}
	for row in csv.DictReader(open(outliers_path, encoding="utf8"), delimiter="	"):
		if row["metric"] == "n_retrocopy_loci" and row["ensembl_gene"] == TP53 and row["assembly"] in PROBOSCIDEA:
			ranks[row["assembly"]] = int(row["rank"])
	check(ranks.get("HLloxAfr4", 999) <= 20, f"TP53 ranks {ranks.get('HLloxAfr4')} of {len(genes)} genes in the African elephant retrocopy outliers (must be top 20)")
	check(ranks.get("HLeleMax1", 999) <= 50, f"TP53 ranks {ranks.get('HLeleMax1')} in the Asian elephant retrocopy outliers (must be in the top 50 list)")
print("all checks passed")
