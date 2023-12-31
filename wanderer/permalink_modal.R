#!/usr/bin/env R


# generates a modal with a give url (string)
generate_modal <- function(url) {
    modal_template <- sprintf('
<button type="button" class="btn btn-default" data-toggle="modal" data-target="#myModal" data-remote="http://example.com">Share current plot (permanent link)</button>

<div id="myModal" class="modal hide fade">
    <div class="modal-header">
        <button type="button" class="close" data-dismiss="modal" aria-hidden="true"> x </button>
        <h3 id="myModalLabel">Share the current plot (permanent link)</h3>
    </div>
    <div class="modal-body">        
        <p>
          %s
        </p>
    </div>
    <div class="modal-footer">
        <button class="btn" data-dismiss="modal" aria-hidden="true">Close</button>
    </div>
</div>
', url)
    
}

## generates an url with an start, end and chromosome for hg19
generate_genome_browser_link <- function(chromosome, start, end) {
    base <- 'http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&position=%s%%3A%s-%s'
    loc <- sprintf(base, chromosome, start, end)

    template <- sprintf('<a class="btn" target = "_blank" href="%s">Go to Genome Browser</a>', loc)
    ## template <- sprintf('<button type="button" class="btn" href="%s">Go to Genome Browser</button>',
                        ## loc)                          
}


generate_modal <- function(url) {
    modal_template <- sprintf('
<button type="button" class="btn btn-default" data-toggle="modal" data-target="#myModal" data-remote="http://example.com">Share current plot (permanent link)</button>

<div id="myModal" class="modal hide fade">
    <div class="modal-header">
        <button type="button" class="close" data-dismiss="modal" aria-hidden="true"> x </button>
        <h3 id="myModalLabel">Share the current plot (permanent link)</h3>
    </div>
    <div class="modal-body">        
        <p>
          %s
        </p>
        <p>
          <a href="%s" target="_blank">Follow the link (opens a new tab)</a>
        <p>
    </div>
    <div class="modal-footer">
        <button class="btn" data-dismiss="modal" aria-hidden="true">Close</button>
    </div>
</div>
', url, url)
    
}

pop_modal_plot <- function(url) {
    modal_template <- '
<button type="button" class="btn btn-default" data-toggle="modal" data-target="#modal_plot" data-remote="http://example.com">Non proportional plot</button>

<div id="modal_plot" class="modal hide fade">
    <div class="modal-header">
        <button type="button" class="close" data-dismiss="modal" aria-hidden="true"> x </button>
        <h3 id="myModalLabel">Test</h3>
    </div>
    <div class="modal-body">        
      <img src="test.png" id="imagepreview" style="width: 600px;" >
    </div>


    <div class="modal-footer">
        <button class="btn" data-dismiss="modal" aria-hidden="true">Close</button>
    </div>
</div>
'
}


## selected regulome explorer datasets May 2015
regulome_explorer_datasets <- function() {
    dat <- list('acc'  = 'acc_2015_03_31',
                'blca' = 'blca_20may13_test',
                'brca' = 'brca_manuscript_rerun_nov12d_pw',
                'chol' = 'chol_20150515_private',
                'coad' = 'coad_03feb13_seq_tumor_only',
                'gbm'  = 'gbm_2013_pub_tumor_only',
                'hnsc' = 'hnsc_03feb13_seq_tumor_only',
                'kirc' = 'kirc_01oct12_A_pw',
                'lgg'  = 'lgg_04oct13_seq',
                'luad' = 'luad_03feb13_seq_tumor_only',
                'lusc' = 'lusc_03feb13_seq_tumor_only',
                'ov'   = 'ov_03feb13_ary_tumor_only',
                'skcm' = 'skcm_01apr14_331_all',
                'stad' = 'stad_23jan14_seq_tumor_only',
                'thca' = 'thca_18oct14_TP',
                'ucec' = 'ucec_28jun13b_seq_tumor_only',
                'uvm'  = 'uvm_20150515_private')
    return(dat)
}
 
## generates an url with an start, end and chromosome for hg19
generate_regulome_explorer_link <- function(dataset, gene) {
    dataset <- regulome_explorer_datasets()[[dataset]]
    
    ## 'http://explorer.cancerregulome.org/all_pairs/?dataset=chol_20150515_private&t_type=*&t_label=TP53'
    base <- 'http://explorer.cancerregulome.org/all_pairs/?dataset=%s&t_type=*&t_label=%s'
    
    loc <- sprintf(base, dataset, gene)

    template <- sprintf('<a class="btn" target = "_blank" href="%s">Go to Regulome explorer</a>', loc)
    ## template <- sprintf('<button type="button" class="btn" href="%s">Go to Genome Browser</button>',
                        ## loc)

    return(template)
}


## selected cbioportal datasets May 2015
cbioportal_datasets <- function(x){
    dat <- list(acc = c(rep('acc_tcga', 4), 'acc_tcga_cnaseq'),
                blca = c(rep('blca_tcga_pub', 4), 'acc_tcga_cnaseq'),
                brca = c(rep('brca_tcga', 4), 'brca_tcga_all'),
                cesc = c(rep('cesc_tcga', 4), 'cesc_tcga_cnaseq'),                         
                chol = c(rep('chol_nccs_2013', 4), 'chol_nccs_2013_cnaseq'),
                coad = c(rep('coadread_tcga', 4), 'coadread_tcga_all'),
                esca = c(rep('esca_tcga', 4), 'esca_tcga_cnaseq'),
                gbm = c(rep('gbm_tcga', 4), 'gbm_tcga_all'),                
                hnsc = c(rep('hnsc_tcga', 4), 'hnsc_tcga_cnaseq'),
                kich = c(rep('kich_tcga_pub', 4), 'kich_tcga_pub_cnaseq'),
                kirc = c(rep('kirc_tcga_pub', 4), 'kirc_tcga_pub_cnaseq'),
                kirp = c(rep('kirp_tcga', 4), 'kirp_tcga_cnaseq'),
                laml = c(rep('laml_tcga_pub', 4), 'laml_tcga_pub_cnaseq'),
                lgg = c(rep('lgg_tcga', 4), 'lgg_tcga_cnaseq'),
                lihc = c(rep('lihc_tcga', 4), 'lihc_tcga_cnaseq'),
                luad = c(rep('luad_tcga_pub', 4), 'luad_tcga_all'),                
                lusc = c(rep('lusc_tcga', 4), 'lusc_tcga_cnaseq'),
                ov = c(rep('ov_tcga_pub', 4), 'ov_tcga_pub_cnaseq'),
                paad = c(rep('paad_tcga', 4), 'paad_tcga_cnaseq'),
                pcpg = c(rep('pcpg_tcga', 4), 'pcpg_tcga_cnaseq'),
                prad = c(rep('prad_tcga', 4), 'prad_tcga_cnaseq'), 
                read = c(rep('coadread_tcga', 4), 'coadread_tcga_cnaseq'),
                sarc = c(rep('sarc_tcga', 4), 'sarc_tcga_cnaseq'),
                skcm = c(rep('skcm_tcga', 4), 'skcm_tcga_cnaseq'),
                stad = c(rep('stad_tcga', 4), 'stad_tcga_cnaseq'),
                thca = c(rep('thca_tcga', 4), 'thca_tcga_cnaseq'),
                ucec = c(rep('ucec_tcga', 4), 'ucec_tcga_cnaseq'),
                ucs = c(rep('ucs_tcga', 4), 'ucsc_tcga_cnaseq'))
    return(dat)
}



## http://www.cbioportal.org/index.do?cancer_study_list=acc_tcga&cancer_study_id=acc_tcga&genetic_profile_ids_PROFILE_MUTATION_EXTENDED=acc_tcga_mutations&genetic_profile_ids_PROFILE_COPY_NUMBER_ALTERATION=acc_tcga_gistic&Z_SCORE_THRESHOLD=2.0&RPPA_SCORE_THRESHOLD=2.0&data_priority=0&case_set_id=acc_tcga_cnaseq&case_ids=&gene_set_choice=user-defined-list&gene_list=ACTA1&clinical_param_selection=null&tab_index=tab_visualize&Action=Submit

## http://www.cbioportal.org/index.do?cancer_study_list=coadread_tcga_pub&cancer_study_id=coadread_tcga_pub&genetic_profile_ids_PROFILE_MUTATION_EXTENDED=coadread_tcga_pub_mutations&Z_SCORE_THRESHOLD=2.0&data_priority=0&case_set_id=coadread_tcga_pub_sequenced&case_ids=&gene_set_choice=user-defined-list&gene_list=HDAC9&clinical_param_selection=null&tab_index=tab_visualize&Action=Submit
generate_cbioportal_link <- function(dataset, gene) {
    dataset <- cbioportal_datasets()[[dataset]]

        ## 'http://explorer.cancerregulome.org/all_pairs/?dataset=chol_20150515_private&t_type=*&t_label=TP53'
    base <- 'http://www.cbioportal.org/index.do?cancer_study_list=%s&cancer_study_id=%s&genetic_profile_ids_PROFILE_MUTATION_EXTENDED=%s_mutations&genetic_profile_ids_PROFILE_COPY_NUMBER_ALTERATION=%s_gistic&Z_SCORE_THRESHOLD=2.0&RPPA_SCORE_THRESHOLD=2.0&data_priority=0&case_set_id=%s&case_ids=&gene_set_choice=user-defined-list&gene_list=%s&clinical_param_selection=null&tab_index=tab_visualize&Action=Submit'
    
    loc <- sprintf(base, dataset[1], dataset[2], dataset[3], dataset[4], dataset[5], gene)

    template <- sprintf('<a class="btn" target = "_blank" href="%s">Go to cBioPortal</a>', loc)
    ## template <- sprintf('<button type="button" class="btn" href="%s">Go to Genome Browser</button>',
                        ## loc)                         

    return(template)
}


## http://www.cbioportal.org/index.do?cancer_study_list=cesc_tcga&cancer_study_id=cesc_tcga&genetic_profile_ids_PROFILE_MUTATION_EXTENDED=cesc_tcga_mutations&genetic_profile_ids_PROFILE_COPY_NUMBER_ALTERATION=cesc_tcga_gistic&Z_SCORE_THRESHOLD=2.0&RPPA_SCORE_THRESHOLD=2.0&data_priority=0&case_set_id=cesc_tcga_cnaseq&case_ids=&gene_set_choice=user-defined-list&gene_list=ACTA1&clinical_param_selection=null&tab_index=tab_visualize&Action=Submit
