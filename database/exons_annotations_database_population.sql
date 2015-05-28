-- Populating the exons annotations that have a proper rnaseqgeneid


--backup
create table annotations.exons_annot_bak AS (
       select * from annotations.exons_annot);

drop table annotations.exons_annot;

CREATE TABLE annotations.exons_annot(
id               varchar, 
exon             varchar,     
emsemblgeneid    varchar, 
emsembltransid   varchar, 
chr              varchar,     
exon_start       integer,               
exon_end         integer,               
strand           varchar, 
genestart        integer,               
geneend          integer,               
genebiotype      varchar,                  
genename         varchar,
RNAseqGeneID     varchar)
;

\copy annotations.exons_annot from /imppc/labs/maplab/share/anna2izaskun/db_region_profile_data/exons_annotations_new.tab  WITH null AS 'NA' DELIMITER E'\t'

comment on table annotations.exons_annot is 'Data from Anna Diez';

-- removing the backup
drop table annotations.exons_annot_bak;
