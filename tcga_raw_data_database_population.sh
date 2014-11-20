#!/bin/bash
#
# Fetches TCGA-Assembler-based downloaded data and populates the postgres with it
#
# Table declarations by Anna are badly formatted, on the fly building of proper ones
# the problems are the schema/name tuples
# and even the column names, that containg dots and '-' signs
#
# Izaskun Mallona
# 13 nov 2014

# database stuff

#HOST=172.19.5.24
HOST=overlook
PORT=5432
USER=imallona
DB=regional_profile
PASS=""

# data origins stuff
WDIR=/imppc/labs/maplab/share/anna2izaskun/db_region_profile_data
DAT=(BLCA BRCA CESC COAD GBM  HNSC KIRC KIRP LIHC LUAD LUSC PAAD PCPG PRAD READ SARC SKCM THCA UCEC)

PLATFORM=(450KMeth RNAseq)
CTYPE=(Normal Tumor)


# BLCA_RNAseq_Normal.csv
# BLCA_RNAseq_Normal_createTable.txt

cd $WDIR

for name in ${DAT[@]}
do
    for platform in ${PLATFORM[@]}
    do
        for ctype in ${CTYPE[@]}
        do
            curr_data_fn="$name"/"$type"/"$name"_"$platform"_"$ctype".csv
            curr_sql_fn="$WDIR"/"$name"/"$platform"/"$ctype"/"$name"_"$platform"_"$ctype"_createTable.txt

            # fixing the table declarations start

            awk '{if (NR > 1) print}' "$curr_sql_fn" > /tmp/f

            table="illuminahiseq_rnaseqv2"
            pk_field="exon"
            if [ "$platform" == "450KMeth" ]
            then
                table="humanmethylation450"
                pk_field="probe"
            fi

            schema_table="$name"_"$ctype"."$table"
            header="CREATE TABLE ""$schema_table""("

            echo $header > /tmp/f2
            cat /tmp/f2 /tmp/f > /tmp/stmt
            
            # removing the odd characters
            sed 's/\./_/g' /tmp/stmt > /tmp/stmt2
            sed 's/-/_/g' /tmp/stmt2 > /tmp/stmt
            
            # fixing the table declarations start

            # PGPASSWORD="$PASS" psql -h $HOST -p $PORT -U $USER \
            # -d $DB -t -A -F"," -f "$curr_sql_fn"

            PGPASSWORD="$PASS" psql -h $HOST -p $PORT -U $USER \
            -d $DB -t -A -F"," -f /tmp/stmt

            rm /tmp/stmt /tmp/stmt2 /tmp/f /tmp/f2
            # PGPASSWORD="$PASS" psql -h $HOST -p $PORT -U $USER \
            # -d $DB -t -A -F"," -c "\copy thca_tumor.humanmethylation450 FROM '/imppc/labs/maplab/share/anna2izaskun/db_region_profile_data/thca_tumor/thca_humanmethylation450_to_SQL.csv' WITH NULL AS 'NA' DELIMITER ','"

            # primary key generation start

            echo "CREATE UNIQUE INDEX ""$table"_idx ON "$schema_table"" (""$pk_field"");" > /tmp/stmt_pk
            echo "ALTER TABLE ""$schema_table" >> /tmp/stmt_pk
            echo " ADD CONSTRAINT ""$table""_pkey PRIMARY KEY USING INDEX ""$table""_idx;" >> /tmp/stmt_pk

            # primary key generation end

            rm /tmp/stmt_pk
            
        done
    done
done
