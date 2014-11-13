#!/bin/bash
#
# Fetches TCGA-Assembler-based downloaded data and populates the postgres with it
#
# Izaskun Mallona
# 13 nov 2014

# database stuff

HOST=172.19.5.24
# HOST=overlook
PORT=5432
USER=imallona
DB=meth_correlations
PASS=""

# data origins stuff
WDIR=/imppc/labs/maplab/share/anna2izaskun/db_region_profile_data
DAT=(BLCA, BRCA, CESC, COAD, GBM,  HNSC, KIRC, KIRP, LIHC, LUAD, LUSC, PAAD, PCPG, PRAD, READ, SARC, SKCM, THCA, UCEC)

PLATFORM=(450KMeth,RNAseq)
CTYPE=(Normal,Tumor)


# BLCA_RNAseq_Normal.csv
# BLCA_RNAseq_Normal_createTable.txt

cd $WDIR

for name in ${DAT[@]}
do
    for platform in {PLATFORM[@]}
    do
        for ctype in {CTYPE[@]}
        do
            curr_data_fn="$name"/"$type"/"$name"_"$platform"_"$ctype".csv
            curr_sql_fn="$name"/"$type"/"$name"_"$platform"_"$ctype"_createTable.txt
          
            PGPASSWORD="$PASS" psql -h $HOST -p $PORT -U $USER \
            -d $DB -t -A -F"," -f "$curr_sql_fn"

            # PGPASSWORD="$PASS" psql -h $HOST -p $PORT -U $USER \
            # -d $DB -t -A -F"," -c "\copy thca_tumor.humanmethylation450 FROM '/imppc/labs/maplab/share/anna2izaskun/db_region_profile_data/thca_tumor/thca_humanmethylation450_to_SQL.csv' WITH NULL AS 'NA' DELIMITER ','"
            
        done
    done
done
