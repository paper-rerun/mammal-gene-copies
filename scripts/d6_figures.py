# Static figures for the dataset (PNG at 2x + SVG, light and dark variants). Every number comes from out/.
# Usage: d6_figures.py <out_dir> <anage_data.txt> <fig_dir>
import csv, gzip, math, os, sys
from collections import defaultdict
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

out_dir, anage_path, fig_dir = sys.argv[1:4]
os.makedirs(fig_dir, exist_ok=True)

# Reference palette (dataviz skill, references/palette.md), light and dark steps
THEMES = {
	"light": {"surface": "#fcfcfb", "ink": "#0b0b0b", "ink2": "#52514e", "muted": "#898781", "grid": "#e1e0d9", "axis": "#c3c2b7",
		"s1": "#2a78d6", "s2": "#eb6834", "s3": "#1baf7a", "seq": ["#fcfcfb", "#cde2fb", "#86b6ef", "#3987e5", "#1c5cab", "#0d366b"],
		"div_neg": "#2a78d6", "div_mid": "#f0efec", "div_pos": "#e34948", "background_dots": "#c3c2b7"},
	"dark": {"surface": "#1a1a19", "ink": "#f0efec", "ink2": "#c3c2b7", "muted": "#898781", "grid": "#2c2c2a", "axis": "#383835",
		"s1": "#3987e5", "s2": "#d95926", "s3": "#199e70", "seq": ["#1a1a19", "#184f95", "#256abf", "#3987e5", "#86b6ef", "#cde2fb"],
		"div_neg": "#3987e5", "div_mid": "#383835", "div_pos": "#e66767", "background_dots": "#52514e"},
}
TP53 = "ENSG00000141510"

def style(theme):
	t = THEMES[theme]
	plt.rcParams.update({"figure.facecolor": t["surface"], "axes.facecolor": t["surface"], "savefig.facecolor": t["surface"],
		"text.color": t["ink"], "axes.labelcolor": t["ink2"], "xtick.color": t["muted"], "ytick.color": t["muted"],
		"axes.edgecolor": t["axis"], "axes.grid": True, "grid.color": t["grid"], "grid.linewidth": 0.6, "grid.linestyle": "-",
		"axes.spines.top": False, "axes.spines.right": False, "font.family": ["Segoe UI", "DejaVu Sans"], "font.size": 9,
		"axes.titlesize": 10, "axes.titleweight": "bold", "axes.titlelocation": "left", "axes.axisbelow": True})
	return t

def save(fig, name, theme):
	for extension, kwargs in (("png", {"dpi": 200}), ("svg", {})):
		fig.savefig(os.path.join(fig_dir, f"{name}_{theme}.{extension}"), bbox_inches="tight", **kwargs)
	plt.close(fig)

# ---------- data ----------
species = list(csv.DictReader(open(os.path.join(out_dir, "species.tsv"), encoding="utf8"), delimiter="\t"))
by_assembly = {s["assembly"]: s for s in species}
lineage_of = {r["assembly"]: r["lineage"] for r in csv.DictReader(open(os.path.join(out_dir, "toga_species.tsv"), encoding="utf8"), delimiter="\t")}
genes = {g["ensembl_gene"]: g for g in csv.DictReader(open(os.path.join(out_dir, "genes.tsv"), encoding="utf8"), delimiter="\t")}
LOSS_GENES = {"UCP1": "ENSG00000109424", "F12": "ENSG00000131187", "KLKB1": "ENSG00000164344", "MMP20": "ENSG00000137674"}
wanted = set(LOSS_GENES.values()) | {TP53}
cell = {}
with gzip.open(os.path.join(out_dir, "gene_by_species.tsv.gz"), "rt", encoding="utf8") as f:
	next(f)
	for line in f:
		assembly, gene, status, _, n_loci, n_retro = line.rstrip("\n").split("\t")
		if gene in wanted:
			cell[(assembly, gene)] = (status, n_loci, n_retro)

def clade(assembly):
	lineage = lineage_of[assembly]
	if "Cetacea" in lineage:
		return "Cetacea"
	if "Proboscidea" in lineage:
		return "Proboscidea"
	if "Sirenia" in lineage:
		return "Sirenia"
	if "Hippopotamidae" in lineage:
		return "Hippopotamidae"
	return by_assembly[assembly]["toga_order_folder"]

def short(name):
	parts = name.split()
	return parts[0][0] + ". " + parts[1] if len(parts) > 1 else name

# ---------- figure 1: TP53 retrocopies per assembly ----------
def figure_tp53(theme):
	t = style(theme)
	# Sparse retrocopy scans are left out: their zeros are not comparable (qc_retrocopy_scan_sparse)
	rows = [(int(v[2]), a) for (a, g), v in cell.items() if g == TP53 and v[2] != "NA" and by_assembly[a]["qc_retrocopy_scan_sparse"] == "0"]
	rows.sort()
	values = np.array([r[0] for r in rows])
	labels = {"HLloxAfr4": "African elephant", "HLeleMax1": "Asian elephant", "HLbalMys1": "bowhead whale",
		"HLhetGla3": "naked mole-rat", "myoBra1": "Brandt's bat", "triMan1": "Florida manatee"}
	fig, ax = plt.subplots(figsize=(8, 3.6))
	x = np.arange(len(rows))
	highlight = np.array([r[1] in ("HLloxAfr4", "HLeleMax1") for r in rows])
	ax.scatter(x[~highlight], values[~highlight], s=9, color=t["background_dots"], linewidths=0, label="other assemblies")
	ax.scatter(x[highlight], values[highlight], s=36, color=t["s2"], edgecolors=t["surface"], linewidths=1.5, label="elephants")
	median = float(np.median(values))
	ax.axhline(median, color=t["ink2"], linewidth=1)
	ax.text(2, median + 0.4, f"median {median:g}", color=t["ink2"], fontsize=8, va="bottom")
	for k, (value, assembly) in enumerate(rows):
		if assembly in labels:
			ax.annotate(f"{labels[assembly]} ({value})", (k, value), xytext=(-14, -4) if value > 10 else (-6, 14), textcoords="offset points",
				ha="right", va="center" if value > 10 else "bottom", fontsize=8, color=t["ink"], arrowprops={"arrowstyle": "-", "color": t["muted"], "linewidth": 0.6})
	ax.set_ylim(-1, values.max() + 3)
	ax.set_xlabel(f"assemblies with a non-sparse retrocopy scan, sorted (n = {len(rows)})")
	ax.set_ylabel("TP53 retrocopy loci")
	ax.set_xticks([])
	ax.set_title("Positive control: TP53 retrocopies per assembly (TOGA processed pseudogenes)")
	ax.legend(frameon=False, loc="upper left", fontsize=8)
	save(fig, "fig1_tp53_retrocopies", theme)

# ---------- figure 2: known losses by clade ----------
def figure_losses(theme):
	t = style(theme)
	groups = defaultdict(list)
	for s in species:
		groups[clade(s["assembly"])].append(s["assembly"])
	order = ["Cetacea", "Hippopotamidae", "Ruminantia", "Suina", "Tylopoda", "Perissodactyla", "Carnivora", "Pholidota", "Chiroptera",
		"Eulipotyphla", "Primates", "Dermoptera", "Scandentia", "Rodentia", "Lagomorpha", "Proboscidea", "Sirenia", "Afrotheria",
		"Xenarthra", "Metatheria", "Prototheria"]
	order = [o for o in order if o in groups]
	gene_names = list(LOSS_GENES)
	share = np.zeros((len(order), len(gene_names)))
	text = [[""] * len(gene_names) for _ in order]
	for i, group in enumerate(order):
		for j, name in enumerate(gene_names):
			statuses = [cell[(a, LOSS_GENES[name])][0] for a in groups[group]]
			called = [s for s in statuses if s in ("I", "PI", "L")]
			lost = sum(s == "L" for s in called)
			share[i, j] = lost / len(called) if called else np.nan
			text[i][j] = f"{lost}/{len(called)}"
	from matplotlib.colors import LinearSegmentedColormap
	cmap = LinearSegmentedColormap.from_list("seq", t["seq"])
	fig, ax = plt.subplots(figsize=(5.2, 6.4))
	ax.grid(False)
	image = ax.imshow(share, cmap=cmap, vmin=0, vmax=1, aspect="auto")
	for i in range(len(order)):
		for j in range(len(gene_names)):
			value = share[i, j]
			dark_cell = (value >= 0.6) if theme == "light" else (value >= 0.6)
			color = ("#ffffff" if theme == "light" else "#0b0b0b") if dark_cell else t["ink"]
			ax.text(j, i, text[i][j], ha="center", va="center", fontsize=7.5, color=color)
	ax.set_xticks(range(len(gene_names)), gene_names)
	ax.xaxis.tick_top()
	labels = [f"{g} (n={len(groups[g])})" for g in order]
	ax.set_yticks(range(len(order)), labels)
	for tick in ax.get_yticklabels():
		if tick.get_text().startswith("Cetacea"):
			tick.set_fontweight("bold")
			tick.set_color(t["ink"])
	ax.add_patch(plt.Rectangle((-0.5, -0.5), len(gene_names), 1, fill=False, edgecolor=t["s2"], linewidth=2))
	for spine in ax.spines.values():
		spine.set_visible(False)
	bar = fig.colorbar(image, ax=ax, fraction=0.04, pad=0.03)
	bar.set_label("share of assemblies with the gene lost (L)\namong those called I, PI or L", color=t["ink2"])
	bar.outline.set_visible(False)
	ax.set_title("Known gene losses by clade (cell text: lost / called)", pad=28)
	save(fig, "fig2_known_losses", theme)

# ---------- figure 3: assembly quality ----------
def figure_quality(theme):
	t = style(theme)
	n50 = np.array([float(s["contig_n50"]) for s in species])
	frac = np.array([float(s["frac_multicopy"]) for s in species])
	flag = np.array([s["qc_false_duplication"] == "1" for s in species])
	pilot_only = np.array([s["qc_false_duplication_pilot35"] == "1" and s["qc_false_duplication"] != "1" for s in species])
	passes = ~flag & ~pilot_only
	threshold = 3 * float(np.median(frac))
	fig, ax = plt.subplots(figsize=(6.4, 4))
	ax.scatter(n50[passes], frac[passes] * 100, s=12, color=t["s1"], alpha=0.7, linewidths=0, label=f"passes (n={int(passes.sum())})")
	ax.scatter(n50[pilot_only], frac[pilot_only] * 100, s=22, marker="o", facecolors="none", edgecolors=t["s2"], linewidths=1.2,
		label=f"pilot list only (n={int(pilot_only.sum())})")
	ax.scatter(n50[flag], frac[flag] * 100, s=22, marker="D", color=t["s2"], edgecolors=t["surface"], linewidths=0.8,
		label=f"flagged, qc_false_duplication (n={int(flag.sum())})")
	ax.axhline(threshold * 100, color=t["ink2"], linewidth=1)
	ax.text(n50.max(), threshold * 100 + 0.6, f"flag threshold {threshold * 100:.1f}% (3 x median)", ha="right", fontsize=8, color=t["ink2"])
	ax.set_xscale("log")
	ax.set_xlabel("contig N50 (bp, log scale)")
	ax.set_ylabel("genes with 2+ ortholog loci (%)")
	ax.set_title("Assembly quality: apparent duplication vs contig N50")
	ax.legend(frameon=False, fontsize=8, loc="upper right", bbox_to_anchor=(1, 0.92))
	save(fig, "fig3_assembly_quality", theme)

# ---------- figure 4: lifespan vs body mass ----------
def figure_life_history(theme):
	t = style(theme)
	rows = [s for s in species if s["max_longevity_yrs"] != "NA" and s["adult_weight_g"] != "NA"]
	seen = set()
	unique = []
	for s in rows:
		if s["anage_name"] not in seen:
			seen.add(s["anage_name"])
			unique.append(s)
	mass = np.log10([float(s["adult_weight_g"]) for s in unique])
	life = np.log10([float(s["max_longevity_yrs"]) for s in unique])
	slope, intercept = np.polyfit(mass, life, 1)
	residual = life - (slope * mass + intercept)
	human = None
	for row in csv.DictReader(open(anage_path, encoding="utf8", errors="replace"), delimiter="\t"):
		if row["Genus"] == "Homo" and row["Species"] == "sapiens":
			human = (math.log10(float(row["Adult weight (g)"] or row["Body mass (g)"])), math.log10(float(row["Maximum longevity (yrs)"])))
	from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
	cmap = LinearSegmentedColormap.from_list("div", [t["div_neg"], t["div_mid"], t["div_pos"]])
	limit = float(np.max(np.abs(residual)))
	fig, ax = plt.subplots(figsize=(6.6, 4.6))
	points = ax.scatter(10 ** mass, 10 ** life, c=residual, cmap=cmap, norm=TwoSlopeNorm(0, -limit, limit), s=16, edgecolors=t["axis"], linewidths=0.3)
	grid = np.linspace(mass.min(), mass.max(), 50)
	ax.plot(10 ** grid, 10 ** (slope * grid + intercept), color=t["ink2"], linewidth=1)
	ax.set_xscale("log")
	ax.set_yscale("log")
	labels = {"Balaena mysticetus": "bowhead whale", "Heterocephalus glaber": "naked mole-rat", "Myotis brandtii": "Brandt's bat", "Loxodonta africana": "African elephant"}
	for s, m, l in zip(unique, mass, life):
		if s["anage_name"] in labels:
			offset = (-70, 0) if s["anage_name"] == "Balaena mysticetus" else (6, 6)
			ax.annotate(labels[s["anage_name"]], (10 ** m, 10 ** l), xytext=offset, textcoords="offset points", fontsize=8, color=t["ink"], va="center" if offset[0] < 0 else "baseline")
	if human:
		ax.scatter([10 ** human[0]], [10 ** human[1]], s=40, facecolors="none", edgecolors=t["ink"], linewidths=1.2)
		ax.annotate("human (AnAge; not a TOGA query)", (10 ** human[0], 10 ** human[1]), xytext=(6, -12), textcoords="offset points", fontsize=8, color=t["ink"])
	bar = fig.colorbar(points, ax=ax, fraction=0.04, pad=0.02)
	bar.set_label("residual log10 lifespan (above/below the line)", color=t["ink2"])
	bar.outline.set_visible(False)
	ax.set_xlabel("adult body mass (g, log scale)")
	ax.set_ylabel("maximum lifespan (years, log scale)")
	ax.set_title(f"Lifespan vs body mass, {len(unique)} species in the table (OLS line, no phylogeny)")
	save(fig, "fig4_lifespan_vs_mass", theme)

# ---------- figure 5: per-species outliers ----------
def figure_outliers(theme):
	t = style(theme)
	chosen = [("HLloxAfr4", "African elephant"), ("HLeleMax1", "Asian elephant"), ("HLtriManLat2", "West Indian manatee"),
		("HLbalMys1", "bowhead whale"), ("HLbalMus2", "blue whale"), ("HLphyCat2", "sperm whale"),
		("hetGla2", "naked mole-rat"), ("HLfukDam2", "Damara mole-rat"), ("myoBra1", "Brandt's bat"),
		("HLmyoMyo6", "greater mouse-eared bat"), ("mm39", "house mouse"), ("HLpanTro7", "chimpanzee")]
	available = {s["assembly"] for s in species}
	top = defaultdict(list)
	for row in csv.DictReader(open(os.path.join(out_dir, "outliers_by_assembly.tsv"), encoding="utf8"), delimiter="\t"):
		if row["metric"] == "n_retrocopy_loci" and (int(row["rank"]) <= 5 or (row["ensembl_gene"] == TP53 and row["assembly"] in ("HLloxAfr4", "HLeleMax1"))):
			top[row["assembly"]].append(row)
	chosen = [(a, n) for a, n in chosen if a in available and a in top][:12]
	fig, axes = plt.subplots(4, 3, figsize=(9, 9.5), sharex=True)
	z_max = max(float(r["robust_z"]) for a, _ in chosen for r in top[a])
	for ax, (assembly, name) in zip(axes.flat, chosen):
		rows = top[assembly][::-1]
		y = np.arange(len(rows))
		z = [float(r["robust_z"]) for r in rows]
		colors = [t["s2"] if r["ensembl_gene"] == TP53 else t["s1"] for r in rows]
		ax.hlines(y, 0, z, color=t["axis"], linewidth=1)
		ax.scatter(z, y, color=colors, s=30, zorder=3)
		ax.set_yticks(y, [f"{r['toga_symbol']} ({r['value']})" + (f" rank {r['rank']}" if int(r["rank"]) > 5 else "") for r in rows], fontsize=7.5)
		marks = (" *" if by_assembly[assembly]["qc_false_duplication"] == "1" else "") + (" †" if by_assembly[assembly]["qc_retrocopy_scan_sparse"] == "1" else "")
		ax.set_title(name + marks, fontsize=9)
		ax.set_xscale("symlog", linthresh=10)
		ax.set_xlim(0, z_max * 1.3)
		ax.grid(axis="y", visible=False)
	for ax in axes.flat[len(chosen):]:
		ax.set_visible(False)
	for ax in axes[-1]:
		ax.set_xlabel("robust z vs all other assemblies (symlog)")
	fig.suptitle("Top 5 retrocopy outliers per species; label: gene (retrocopy loci). Elephants also show TP53 (orange) with its genome-wide rank.\n"
		"* false-duplication flag   † sparse retrocopy scan (counts not comparable)", x=0.01, ha="left", fontweight="bold", fontsize=9.5)
	fig.tight_layout()
	save(fig, "fig5_outliers", theme)

for theme in THEMES:
	figure_tp53(theme)
	figure_losses(theme)
	figure_quality(theme)
	figure_life_history(theme)
	figure_outliers(theme)
print("figures written to", fig_dir)
