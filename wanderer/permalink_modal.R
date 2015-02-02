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
