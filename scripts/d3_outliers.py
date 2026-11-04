# Per-lineage outlier ranking: for each assembly, genes where its count sits far above all other assemblies.
# robust z = (x - median_others) / (1.4826 * max(MAD_others, 0.5)); metrics: n_retrocopy_loci, n_ortholog_loci.
# Also reports how many other assemblies reach the same value, and how many of those share the TOGA order folder,
# so a clade-wide expansion (e.g. both elephants) can be told apart from a single-assembly one.
# The baseline ("others") leaves out assemblies whose QC flag makes the metric unreliable: sparse retrocopy scans for
# n_retrocopy_loci, false-duplication assemblies for n_ortholog_loci. Flagged assemblies are still ranked themselves.
# Usage: d3_outliers.py <out_dir> [top_n]
import csv, gzip, os, sys
import numpy as np

out_dir = sys.argv[1]
top_n = int(sys.argv[2]) if len(sys.argv) > 2 else 50

species = list(csv.DictReader(open(os.path.join(out_dir, "species.tsv"), encoding="utf8"), delimiter="\t"))
assemblies = [s["assembly"] for s in species]
column_of = {a: k for k, a in enumerate(assemblies)}
order_of = np.array([s["toga_order_folder"] for s in species])
flag_of = {s["assembly"]: s["qc_false_duplication"] for s in species}
sparse_of = {s["assembly"]: s["qc_retrocopy_scan_sparse"] for s in species}
unreliable = {"n_retrocopy_loci": np.array([s["qc_retrocopy_scan_sparse"] == "1" for s in species]),
	"n_ortholog_loci": np.array([s["qc_false_duplication"] == "1" for s in species])}
genes = [g["ensembl_gene"] for g in csv.DictReader(open(os.path.join(out_dir, "genes.tsv"), encoding="utf8"), delimiter="\t")]
symbol_of = {g["ensembl_gene"]: g["toga_symbol"] for g in csv.DictReader(open(os.path.join(out_dir, "genes.tsv"), encoding="utf8"), delimiter="\t")}
row_of = {g: k for k, g in enumerate(genes)}

loci = np.full((len(genes), len(assemblies)), -1, dtype=np.int32)
retro = np.full((len(genes), len(assemblies)), -1, dtype=np.int32)
with gzip.open(os.path.join(out_dir, "gene_by_species.tsv.gz"), "rt", encoding="utf8") as f:
	next(f)
	for line in f:
		assembly, gene, _, _, n_loci, n_retro = line.rstrip("\n").split("\t")
		r, c = row_of[gene], column_of[assembly]
		if n_loci != "NA":
			loci[r, c] = int(n_loci)
		if n_retro != "NA":
			retro[r, c] = int(n_retro)

with open(os.path.join(out_dir, "outliers_by_assembly.tsv"), "w", encoding="utf8") as out:
	out.write("assembly\tmetric\trank\tensembl_gene\ttoga_symbol\tvalue\tmedian_others\tmad_others\trobust_z\tmax_others\t"
		"n_others_at_or_above\tn_same_order_at_or_above\tqc_false_duplication\tqc_retrocopy_scan_sparse\n")
	for metric, matrix in (("n_retrocopy_loci", retro), ("n_ortholog_loci", loci)):
		has_data = (matrix >= 0).any(axis=0)
		for c, assembly in enumerate(assemblies):
			if not has_data[c]:
				continue
			others = np.array([k for k in range(len(assemblies)) if k != c and has_data[k] and not unreliable[metric][k]])
			other_values = matrix[:, others].astype(float)
			other_values[other_values < 0] = np.nan
			center = np.nanmedian(other_values, axis=1)
			mad = np.nanmedian(np.abs(other_values - center[:, None]), axis=1)
			scale = 1.4826 * np.maximum(mad, 0.5)
			value = matrix[:, c].astype(float)
			value[value < 0] = np.nan
			z = (value - center) / scale
			z[np.isnan(z)] = -np.inf
			ranking = np.argsort(-z)[:top_n]
			for rank, r in enumerate(ranking, start=1):
				if not np.isfinite(z[r]):
					break
				at_or_above = other_values[r] >= value[r]
				same_order = int(np.sum(at_or_above & (order_of[others] == order_of[c])))
				out.write(f"{assembly}\t{metric}\t{rank}\t{genes[r]}\t{symbol_of[genes[r]]}\t{int(value[r])}\t{center[r]:g}\t{mad[r]:g}\t{z[r]:.2f}\t"
					f"{int(np.nanmax(other_values[r]))}\t{int(np.sum(at_or_above))}\t{same_order}\t{flag_of[assembly]}\t{sparse_of[assembly]}\n")
		print(metric, "done", flush=True)
