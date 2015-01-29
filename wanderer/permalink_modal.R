#!/usr/bin/env R


# generates a modal with a give url (string)
generate_modal <- function(url) {
    modal_template <- sprintf('
<button type="button" class="btn" data-toggle="modal" data-target="#myModal" data-remote="http://example.com">Share plots</button>

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
