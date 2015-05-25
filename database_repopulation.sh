#!/bin/bash
#
# Fetches TCGA-Assembler-based downloaded data and populates the postgres with it
#
# Table declarations by Anna are badly formatted, on the fly building of proper ones
# the problems are the schema/name tuples
# and even the column names, that containg dots and '-' signs
#
# For rnaseq data summarized by gene and for rnaseq by exon
#
# This version replaces some schemas that were not illumina hiseq datasets (Anna's issue)
#
# Izaskun Mallona
# 25 may 2015

# database stuff
HOST=overlook
PORT=5432
USER=imallona
DB=regional_profile
PASS=""

# data origins stuff
# WDIR=/imppc/labs/maplab/share/anna2izaskun/db_region_profile_data/
WDIR=/imppc/labs/maplab/share/anna2izaskun/db_region_profile_data/fixed

DAT=(COAD READ SKCM UCEC UCS)

PLATFORM=(RNAseqGene RNAseq)
#  Anna did inconsistent filenaming
# PLATFORMFN=(RNAseq_Gene RNAseq)
PLATFORMFN=RNAseq
CTYPE=(Normal Tumor)


# BLCA_RNAseq_Normal.csv
# BLCA_RNAseq_Normal_createTable.txt

cd $WDIR


# for name in ${DAT[@]}
for name in ${DAT[@]}
do
    for ctype in ${CTYPE[@]}
    do
        for platform in ${PLATFORM[@]}
        do
            # catching the data source inconsistent filenaming
            if [ "$platform" = "RNAseqGene" ]
            then
                PLATFORMFN=RNAseq_Gene
            else
                PLATFORMFN=RNAseq
            fi

            curr_data_fn="$WDIR"/"$name"/"$platform"/"$ctype"/"$name"_"$PLATFORMFN"_"$ctype".csv
            curr_sql_fn="$WDIR"/"$name"/"$platform"/"$ctype"/"$name"_"$PLATFORMFN"_"$ctype"_createTable.txt

            #  checking whether the filename exists
            if [ -f $curr_data_fn ] ; then

                # fixing the table declarations start
                awk '{if (NR > 1) print}' "$curr_sql_fn" > /tmp/f
                # removing the odd characters
                sed 's/\./_/g' /tmp/f > /tmp/f2
                sed 's/-/_/g' /tmp/f2 > /tmp/f

                table="illuminahiseq_rnaseqv2_by_gene"
                pk_field="gene"

                schema_table="$name"_"$ctype"."$table"
                # echo $schema_table >> /tmp/flo

                header="CREATE TABLE ""$schema_table""("

                echo $header > /tmp/f2
                cat /tmp/f2 /tmp/f > /tmp/stmt
                
                
                # fixing the table declarations start

                # schema statement start
                echo "CREATE SCHEMA IF NOT EXISTS ""$name"_"$ctype"";" > /tmp/schema_stmt
                echo "GRANT USAGE ON SCHEMA ""$name"_"$ctype"" TO maplabr;" >> /tmp/schema_stmt
                echo "ALTER DEFAULT PRIVILEGES IN SCHEMA ""$name"_"$ctype" >> /tmp/schema_stmt
                echo " GRANT SELECT ON TABLES TO maplabr;" >> /tmp/schema_stmt

                <<EOF            

            PGPASSWORD="$PASS" psql -h $HOST -p $PORT -U $USER \
            -d $DB -t -A -F"," -f /tmp/schema_stmt
EOF

                cat /tmp/schema_stmt
                # schema statement end

                # PGPASSWORD="$PASS" psql -h $HOST -p $PORT -U $USER \
                # -d $DB -t -A -F"," -f "$curr_sql_fn"

                <<EOF


            PGPASSWORD="$PASS" psql -h $HOST -p $PORT -U $USER \
            -d $DB -t -A -F"," -f /tmp/stmt
EOF

                echo "TRUNCATE $schema_table;"
                
                cat /tmp/stmt
                rm /tmp/stmt  /tmp/f /tmp/f2
                # PGPASSWORD="$PASS" psql -h $HOST -p $PORT -U $USER \
                # -d $DB -t -A -F"," -c "\copy thca_tumor.humanmethylation450 FROM '/imppc/labs/maplab/share/anna2izaskun/db_region_profile_data/thca_tumor/thca_humanmethylation450_to_SQL.csv' WITH NULL AS 'NA' DELIMITER ','"

                # primary key generation start

                echo "\copy ""$schema_table"" from ""$curr_data_fn""  WITH null AS 'NA' DELIMITER ','" \
                    >> /tmp/stmt_pk
                
                echo "CREATE UNIQUE INDEX ""$table"_idx ON "$schema_table"" (""$pk_field"");" >> /tmp/stmt_pk
                echo "ALTER TABLE ""$schema_table" >> /tmp/stmt_pk
                echo " ADD CONSTRAINT ""$table""_pkey PRIMARY KEY USING INDEX ""$table""_idx;" >> /tmp/stmt_pk
                <<EOF

            PGPASSWORD="$PASS" psql -h $HOST -p $PORT -U $USER \
            -d $DB -t -A -F"," -f /tmp/stmt_pk
EOF
                # primary key generation end


                cat /tmp/stmt_pk
                rm /tmp/stmt_pk


                # # this does not work as the rnaseqgeneid element is not unique 
                # echo "ALTER TABLE ""$schema_table" >> /tmp/stmt_fk
                # echo "  ADD CONSTRAINT gene_fk FOREIGN KEY (gene) REFERENCES annotations.exons_annot (RNAseqGeneID);" >> /tmp/stmt_fk

                # cat /tmp/stmt_fk
                # rm /tmp/stmt_fk




            else
                echo "--$curr_data_fn is empty, skipping."
            fi ;
            
        done
    done
done
