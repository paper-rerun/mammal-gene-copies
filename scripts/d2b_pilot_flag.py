# Add the pilot's false-duplication list (35 assemblies; >14.5% of the 890 pilot panel genes with 2+ ortholog loci,
# i.e. 3x the panel median of 4.8%) as an extra column, for continuity with PILOTO.md.
# Usage: d2b_pilot_flag.py <species.tsv> <pilot exclude_high_multicopy.txt>
import csv, sys

species_path, pilot_path = sys.argv[1:3]
pilot = set(line.strip() for line in open(pilot_path, encoding="utf8") if line.strip())
rows = list(csv.DictReader(open(species_path, encoding="utf8"), delimiter="\t"))
fields = list(rows[0].keys())
if "qc_false_duplication_pilot35" not in fields:
	fields.append("qc_false_duplication_pilot35")
with open(species_path, "w", encoding="utf8", newline="") as out:
	writer = csv.DictWriter(out, fieldnames=fields, delimiter="\t", lineterminator="\n")
	writer.writeheader()
	for row in rows:
		row["qc_false_duplication_pilot35"] = int(row["species"] in pilot)
		writer.writerow(row)
print("pilot35 flagged:", sum(row["species"] in pilot for row in rows), "| strict flag:", sum(row["qc_false_duplication"] == "1" for row in rows))
