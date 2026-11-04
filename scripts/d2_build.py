# Build the dataset: genes table, species table (traits + assembly QC) and long gene x species table.
# Usage: d2_build.py <pilot_dl_dir> <dataset_dl_dir> <out_dir>
import csv, gzip, os, re, sys
from collections import defaultdict
from statistics import median

pilot_dl, data_dl, out_dir = sys.argv[1:4]
species_rows = list(csv.DictReader(open(os.path.join(out_dir, "toga_species.tsv"), encoding="utf8"), delimiter="\t"))

# ---------- genes table ----------
symbol_by_gene = {}
for row in csv.DictReader(open(os.path.join(out_dir, "reference_transcripts.tsv"), encoding="utf8"), delimiter="\t"):
	symbol_by_gene.setdefault(row["t_gene"], row["symbol"])
genes = sorted(symbol_by_gene)

hgnc_by_ensembl = {}
ensembl_by_symbol = {}
ensembl_by_previous = {}
ensembl_by_entrez = {}
with open(os.path.join(pilot_dl, "genesets", "hgnc_complete_set.txt"), encoding="utf8") as f:
	for row in csv.DictReader(f, delimiter="\t"):
		ensembl_id = row["ensembl_gene_id"]
		if not ensembl_id:
			continue
		hgnc_by_ensembl[ensembl_id] = row
		ensembl_by_symbol[row["symbol"]] = ensembl_id
		for previous in row["prev_symbol"].strip('"').split("|"):
			if previous:
				ensembl_by_previous.setdefault(previous, ensembl_id)
		if row["entrez_id"]:
			ensembl_by_entrez[row["entrez_id"]] = ensembl_id

def to_ensembl(symbol):
	return ensembl_by_symbol.get(symbol) or ensembl_by_previous.get(symbol)

# DNA repair: Wood, Lange & Wood list
dna_repair = set()
html = open(os.path.join(pilot_dl, "genesets", "wood_dna_repair.html"), encoding="utf8", errors="replace").read()
for row_html in re.findall(r"<tr[^>]*>(.*?)</tr>", html, re.S):
	cells = [re.sub(r"<[^>]+>|&nbsp;", " ", c).strip() for c in re.findall(r"<t[dh][^>]*>(.*?)</t[dh]>", row_html, re.S)]
	if len(cells) >= 4 and re.match(r"^N[MR]_", cells[3]):
		ensembl_id = to_ensembl(re.split(r"[\s(]", cells[0])[0])
		if ensembl_id:
			dna_repair.add(ensembl_id)

# TP53 pathway: Reactome R-HSA-3700989 "Transcriptional Regulation by TP53"
tp53_pathway = set()
for line in open(os.path.join(data_dl, "reactome", "ReactomePathways.gmt"), encoding="utf8"):
	parts = line.rstrip("\n").split("\t")
	if parts[1] == "R-HSA-3700989":
		for symbol in parts[2:]:
			ensembl_id = to_ensembl(symbol)
			if ensembl_id:
				tp53_pathway.add(ensembl_id)

# Tumour suppressors: IntOGen compendium, role LoF only (ambiguous allowed)
roles = defaultdict(set)
with open(os.path.join(pilot_dl, "genesets", "2024-06-18_IntOGen-Drivers", "Compendium_Cancer_Genes.tsv"), encoding="utf8") as f:
	for row in csv.DictReader(f, delimiter="\t"):
		roles[row["SYMBOL"]].add(row["ROLE"])
tumour_suppressor = {to_ensembl(s) for s, r in roles.items() if "LoF" in r and "Act" not in r} - {None}

genage = set()
with open(os.path.join(pilot_dl, "hagr", "genage_human.csv"), encoding="utf8", errors="replace") as f:
	for row in csv.DictReader(f):
		ensembl_id = to_ensembl(row["symbol"])
		if ensembl_id:
			genage.add(ensembl_id)

with open(os.path.join(out_dir, "genes.tsv"), "w", encoding="utf8") as out:
	out.write("ensembl_gene\ttoga_symbol\thgnc_symbol\thgnc_locus_group\tset_dna_repair_wood\tset_tp53_reactome\tset_tumour_suppressor_intogen\tset_genage_human\n")
	for gene in genes:
		hgnc = hgnc_by_ensembl.get(gene, {})
		flags = [int(gene in s) for s in (dna_repair, tp53_pathway, tumour_suppressor, genage)]
		out.write(f"{gene}\t{symbol_by_gene[gene]}\t{hgnc.get('symbol', '')}\t{hgnc.get('locus_group', '')}\t" + "\t".join(map(str, flags)) + "\n")
print("genes", len(genes), "| dna_repair", len(dna_repair & set(genes)), "tp53_reactome", len(tp53_pathway & set(genes)),
	"tumour_suppressor", len(tumour_suppressor & set(genes)), "genage", len(genage & set(genes)), flush=True)

# ---------- long table ----------
STATUS_KEEP = {"I", "PI", "UL", "L", "M", "PM", "PG"}
qc = {}
long_path = os.path.join(out_dir, "gene_by_species.tsv.gz")
n_rows = 0
with gzip.open(long_path, "wt", encoding="utf8", compresslevel=6) as long_out:
	long_out.write("assembly\tensembl_gene\ttoga_status\torthology_class\tn_ortholog_loci\tn_retrocopy_loci\n")
	for row in species_rows:
		name = os.path.basename(row["toga_dir"])
		folder = os.path.join(pilot_dl, "toga", name)
		classes = defaultdict(set)
		q_genes = defaultdict(set)
		with gzip.open(os.path.join(folder, "orthologsClassification.tsv.gz"), "rt", encoding="utf8") as f:
			next(f)
			for line in f:
				t_gene, _, q_gene, _, orthology_class = line.rstrip("\n").split("\t")
				classes[t_gene].add(orthology_class)
				if q_gene != "None":
					q_genes[t_gene].add(q_gene)
		status = {}
		with gzip.open(os.path.join(folder, "loss_summ_data.tsv.gz"), "rt", encoding="utf8") as f:
			for line in f:
				if line.startswith("GENE\t"):
					_, gene, value = line.rstrip("\n").split("\t")
					status[gene] = value
		loci_path = os.path.join(out_dir, "pseudo_loci", name + ".loci.tsv.gz")
		retro = None
		if os.path.exists(loci_path):
			retro = defaultdict(int)
			with gzip.open(loci_path, "rt", encoding="utf8") as f:
				next(f)
				for line in f:
					retro[line.split("\t", 1)[0]] += 1
		n_genes = n_multicopy = n_with_ortholog = n_lost = n_called = n_missing = 0
		for gene in genes:
			gene_status = status.get(gene, "NA")
			n_loci = len(q_genes[gene]) if gene in classes else "NA"
			class_text = "|".join(sorted(classes[gene])) if gene in classes else "NA"
			retro_text = "NA" if retro is None else retro.get(gene, 0)
			long_out.write(f"{row['assembly']}\t{gene}\t{gene_status}\t{class_text}\t{n_loci}\t{retro_text}\n")
			n_rows += 1
			n_genes += 1
			if n_loci != "NA" and n_loci >= 1:
				n_with_ortholog += 1
				n_multicopy += n_loci >= 2
			if gene_status in ("I", "PI", "L"):
				n_called += 1
				n_lost += gene_status == "L"
			n_missing += gene_status in ("M", "PM")
		qc[row["assembly"]] = {
			"frac_multicopy": n_multicopy / max(n_with_ortholog, 1),
			"frac_lost": n_lost / max(n_called, 1),
			"frac_missing": n_missing / n_genes,
			"has_retrocopy_scan": int(retro is not None),
			"n_retrocopy_loci_total": sum(retro.values()) if retro is not None else "NA",
		}
		print(row["assembly"], "done", flush=True)
print("long rows", n_rows, flush=True)

# ---------- species table ----------
anage = {}
with open(os.path.join(pilot_dl, "hagr", "anage_data.txt"), encoding="utf8", errors="replace") as f:
	for row in csv.DictReader(f, delimiter="\t"):
		if row["Class"] == "Mammalia":
			anage[row["Genus"] + " " + row["Species"]] = row

overview = {}
with open(os.path.join(pilot_dl, "overview.table.tsv"), encoding="utf8") as f:
	for row in csv.DictReader(f, delimiter="\t"):
		overview[row["Assembly name"].strip()] = row

# Upham et al. 2019 tree tips by "Genus species"
tree_text = open(os.path.join(pilot_dl, "tree", "MamPhy_DNAonly_4098sp_NDexp_MCC.tre"), encoding="utf8").read()
tips = {}
for tip in re.findall(r"([A-Z][a-z]+_[a-z]+_[A-Z]+_[A-Z]+)", tree_text):
	parts = tip.split("_")
	tips.setdefault(parts[0] + " " + parts[1], tip)

def number(text):
	try:
		return float(text)
	except (TypeError, ValueError):
		return None

median_multicopy = median(q["frac_multicopy"] for q in qc.values())
threshold = 3 * median_multicopy
# Retrocopy scan regime: ~100 assemblies have a TOGA processed-pseudogene file with very few loci
# (log10 histogram has a separate cluster at 10^2.75-10^3.75); below 10,000 merged loci the counts are not comparable
SPARSE_RETRO_LOCI = 10000
with open(os.path.join(out_dir, "species.tsv"), "w", encoding="utf8") as out:
	out.write("assembly\tspecies\tcommon_name\ttoga_order_folder\tncbi_accession_or_source\tcontig_n50\tscaffold_n50\t"
		"anage_name\tmax_longevity_yrs\tadult_weight_g\tanage_data_quality\tanage_sample_size\tupham2019_tip\t"
		"has_retrocopy_scan\tn_retrocopy_loci_total\tqc_retrocopy_scan_sparse\tfrac_multicopy\tfrac_lost\tfrac_missing\tqc_false_duplication\n")
	for row in species_rows:
		assembly = row["assembly"]
		info = overview.get(assembly, {})
		binomial = " ".join(row["species"].split()[:2])
		record = anage.get(row["species"]) or anage.get(binomial)
		longevity = number(record["Maximum longevity (yrs)"]) if record else None
		weight = (number(record["Adult weight (g)"]) or number(record["Body mass (g)"])) if record else None
		anage_name = (record["Genus"] + " " + record["Species"]) if record else ""
		tip = tips.get(anage_name) or tips.get(binomial) or ""
		q = qc[assembly]
		flag = int(q["frac_multicopy"] > threshold)
		values = [assembly, row["species"], row["common_name"], row["toga_dir"].split("/")[0], info.get("NCBI accession / source", "").strip(),
			row["contig_n50"], info.get("scaffold N50 (bp)", "").strip(), anage_name,
			"NA" if longevity is None else longevity, "NA" if weight is None else weight,
			record["Data quality"] if record else "", record["Sample size"] if record else "", tip,
			q["has_retrocopy_scan"], q["n_retrocopy_loci_total"],
			"NA" if q["n_retrocopy_loci_total"] == "NA" else int(q["n_retrocopy_loci_total"] < SPARSE_RETRO_LOCI), round(q["frac_multicopy"], 5), round(q["frac_lost"], 5), round(q["frac_missing"], 5), flag]
		out.write("\t".join(str(v) for v in values) + "\n")
print("qc: sparse retrocopy scans (< 10000 loci):", sum(1 for q in qc.values() if q["n_retrocopy_loci_total"] != "NA" and q["n_retrocopy_loci_total"] < SPARSE_RETRO_LOCI))
print(f"qc: median frac_multicopy {median_multicopy:.5f} threshold (3x median) {threshold:.5f} flagged {sum(q['frac_multicopy'] > threshold for q in qc.values())}")
