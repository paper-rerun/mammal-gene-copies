# Pick one TOGA (hg38 reference) assembly per species: the one with the highest contig N50
# Usage: s2_select_assemblies.py <overview.table.tsv> <toga_dirs.txt> <out.tsv>
import csv, sys

overview_path, dirs_path, out_path = sys.argv[1:4]
n50_by_assembly = {}
species_by_assembly = {}
common_by_assembly = {}
lineage_by_assembly = {}
with open(overview_path, encoding="utf8") as f:
	for row in csv.DictReader(f, delimiter="\t"):
		assembly = row["Assembly name"].strip()
		try:
			n50_by_assembly[assembly] = int(row["contig N50 (bp)"])
		except ValueError:
			n50_by_assembly[assembly] = 0
		species_by_assembly[assembly] = row["Species"].strip()
		common_by_assembly[assembly] = row["Common name"].strip()
		lineage_by_assembly[assembly] = row["Taxonomic Lineage"].strip()

best_by_species = {}
with open(dirs_path, encoding="utf8") as f:
	for line in f:
		path = line.strip().rstrip("/")
		if not path:
			continue
		assembly = path.split("__")[-1]
		if assembly not in species_by_assembly:
			print("not in overview:", path, file=sys.stderr)
			continue
		species = species_by_assembly[assembly]
		# Skip human itself if present
		if species == "Homo sapiens":
			continue
		current = best_by_species.get(species)
		if current is None or n50_by_assembly[assembly] > n50_by_assembly[current[1]]:
			best_by_species[species] = (path, assembly)

with open(out_path, "w", encoding="utf8") as f:
	f.write("species\tcommon_name\tassembly\tcontig_n50\ttoga_dir\tlineage\n")
	for species in sorted(best_by_species):
		path, assembly = best_by_species[species]
		f.write(f"{species}\t{common_by_assembly[assembly]}\t{assembly}\t{n50_by_assembly[assembly]}\t{path}\t{lineage_by_assembly[assembly]}\n")
print(len(best_by_species), "species")
