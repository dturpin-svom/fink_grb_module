#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  4 17:36:25 2021

@author: Damien Turpin : damien.turpin@cea.fr
"""

import numpy as np

def cut_grb_mod(pdf,grb_config):
        
    # Build the filter on the source classification
    src_type_filt = grb_config['src_type_filt']
    cuts_src_type = []
    for cut in src_type_filt:
        cuts_src_type.append(pdf['v:classification'] == cut)
    if len(cuts_src_type) >= 2:
        k = 0
        while k < len(cuts_src_type):
            if k < 2:
                srctype_mask_filter = np.logical_or(cuts_src_type[k],cuts_src_type[k+1])
                k = k+2
            else:
                srctype_mask_filter = np.logical_or(srctype_mask_filter,cuts_src_type[k])
                k = k+1
    else:
        srctype_mask_filter = cuts_src_type[0]
        
    # Build the filter on the Real Bogus classification  
    rb_mask_filter = pdf['i:drb'] >= 0.9
    
    # Build the filter on the serendipituous association between the GRB and the ZTF candidates
    sigma_grb_ass = grb_config['sigma_grb_ass']
    pser_mask_filter = pdf['v:grbSigmaAss'] >= sigma_grb_ass

    pdf = pdf.loc[rb_mask_filter & srctype_mask_filter & pser_mask_filter]
    return pdf