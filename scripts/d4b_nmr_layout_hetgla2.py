# Layout of naked mole-rat retrocopy loci in the hetGla2 assembly (no sequence: hetGla2 genome not downloaded)
# Usage: d4b_nmr_layout_hetgla2.py <out_dir>
import gzip, os, sys
from collections import defaultdict

out_dir = sys.argv[1]
GENES = {"ENSG00000140262": "TCF12", "ENSG00000104320": "NBN", "ENSG00000196628": "TCF4", "ENSG00000111640": "GAPDH"}
loci = defaultdict(list)
path = os.path.join(out_dir, "pseudo_loci", "Heterocephalus_glaber__naked_mole-rat__hetGla2.loci.tsv.gz")
with gzip.open(path, "rt", encoding="utf8") as f:
	next(f)
	for line in f:
		gene, chrom, start, end, _, _ = line.rstrip("\n").split("\t")
		if gene in GENES:
			loci[GENES[gene]].append((chrom, int(start), int(end)))
print("hetGla2 (layout only; no sequence downloaded)")
for name, items in loci.items():
	items.sort()
	per_scaffold = defaultdict(int)
	for chrom, _, _ in items:
		per_scaffold[chrom] += 1
	gaps = [b[1] - a[2] for a, b in zip(items, items[1:]) if a[0] == b[0]]
	lengths = sorted(end - start for _, start, end in items)
	print(f"{name}: {len(items)} loci on {len(per_scaffold)} scaffolds (max per scaffold {max(per_scaffold.values())}); "
		f"length median {lengths[len(lengths) // 2]} min {lengths[0]} max {lengths[-1]}; same-scaffold gaps < 50 kb {sum(g < 50000 for g in gaps)} of {len(gaps)}")
