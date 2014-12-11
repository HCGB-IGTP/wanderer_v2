#!/bin/bash


BASE=/imppc/labs/maplab/share/anna2izaskun/db_region_profile_data/


awk '{OFS=FS="\t"; print $1,"chr"$2,$3,$4,$5, toupper($6)}' "$BASE"/gene_annotations.tab > "$BASE"/gene_annotations_modified.tab


echo "CREATE TABLE annotations.reduced_genome(" > /tmp/stmt
echo "  EmsemblGeneID VARCHAR(32) PRIMARY KEY," >> /tmp/stmt
echo "  chrom VARCHAR(32), " >> /tmp/stmt
echo "  gene_start INTEGER, " >> /tmp/stmt
echo "  gene_end INTEGER, " >> /tmp/stmt
echo "  strand VARCHAR(2)," >> /tmp/stmt
echo "  GeneName VARCHAR(32))" >> /tmp/stmt
echo ";" >> /tmp/stmt

echo "CREATE UNIQUE INDEX reduced_genome_gene_symbol_idx " >> /tmp/stmt
echo "  ON annotations.reduced_genome  (GeneName, EmsemblGeneID)" >> /tmp/stmt
echo ";" >> /tmp/stmt

echo "\copy annotations.reduced_genome from ""$BASE/gene_annotations_modified.tab""  WITH null AS 'NA' DELIMITER E'\t'"  >> /tmp/stmt


# produces 

# CREATE TABLE annotations.reduced_genome(
#   EmsemblGeneID VARCHAR(32) PRIMARY KEY,
#   chrom VARCHAR(32), 
#   gene_start INTEGER, 
#   gene_end INTEGER, 
#   strand VARCHAR(2),
#   GeneName VARCHAR(32))
# ;
# CREATE UNIQUE INDEX reduced_genome_gene_symbol_idx 
#   ON annotations.reduced_genome  (GeneName, EmsemblGeneID)
# ;
# \copy annotations.reduced_genome from /imppc/labs/maplab/share/anna2izaskun/db_region_profile_data//gene_annotations_modified.tab  WITH null AS 'NA' DELIMITER E'\t'

# comment on table annotations.reduced_genome is 'Data from Anna Diez'
