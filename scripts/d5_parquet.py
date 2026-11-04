# Write a Parquet copy of the long table (same columns, typed; NA -> null)
# Usage: d5_parquet.py <out_dir>
import os, sys
import pyarrow as pa
import pyarrow.csv as pacsv
import pyarrow.parquet as pq

out_dir = sys.argv[1]
source = os.path.join(out_dir, "gene_by_species.tsv.gz")
types = {"assembly": pa.string(), "ensembl_gene": pa.string(), "toga_status": pa.string(), "orthology_class": pa.string(),
	"n_ortholog_loci": pa.int16(), "n_retrocopy_loci": pa.int32()}
table = pacsv.read_csv(source,
	parse_options=pacsv.ParseOptions(delimiter="\t"),
	convert_options=pacsv.ConvertOptions(column_types=types, null_values=["NA"], strings_can_be_null=True))
for name in ("assembly", "ensembl_gene", "toga_status", "orthology_class"):
	index = table.schema.get_field_index(name)
	table = table.set_column(index, name, table.column(name).dictionary_encode())
pq.write_table(table, os.path.join(out_dir, "gene_by_species.parquet"), compression="zstd")
print("rows", table.num_rows)
