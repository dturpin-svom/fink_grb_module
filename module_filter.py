#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  4 17:36:25 2021

@author: Damien Turpin : damien.turpin@cea.fr
"""

import numpy as np
from module_utils import get_grb_config, get_pdf_measure
from astro_calc import dc_mag
import pandas as pd
from scipy import special
from math import sqrt

def det_history_filter(list_pdf_obj,grb_date):
    hist_filter = []
    for pdf in list_pdf_obj:
        tstart_filter = pdf['i:jd'] < grb_date.jd
        tend_filter = pdf['i:jd'] > grb_date.jd+90
        det_tag = pdf['i:sigmapsf'] > 0
        
        hist_filter.append(((tstart_filter | tend_filter) & det_tag).any())
        
    hist_filter = [not elem for elem in hist_filter]
    return hist_filter

def src_type_filter(pdfs,grb_config_path):
    
    grb_config = get_grb_config(grb_config_path)
    
    srctypes = grb_config['src_type_filt']
    src_type_filter = pdfs['v:classification'].\
        apply(lambda x: any([k in x for k in srctypes]))
    
    return src_type_filter
    

def cut_grb_mod(pdfs,grb_config_path):
    
    grb_config = get_grb_config(grb_config_path)
        
    # Build the filter mask on the rb score
    rb_threshold = grb_config['rb_threshold']
    
    # get the source type filter mask
    # srctype_filtermask = src_type_filter(pdfs,grb_config_path)
        
    # Build the filter on the Real Bogus classification  
    rb_mask_filter = pdfs['i:drb'] >= rb_threshold
    
    # Build the filter on the serendipituous association between the GRB and the ZTF candidates
    
    pser_mask_filter = pdfs['v:Pser_grb'] >= special.erf(grb_config['sigma_grb_ass']/sqrt(2))

    pdf_filtered = pdfs.loc[rb_mask_filter & pser_mask_filter]
    
    return pdf_filtered

    
def ft_finder(list_pdf_objs,grb_date):
    """

    Search for fast-fading transients in the positive GRB/ZTF crossmatched list
    ----------
    list_pdf_objs : panda data frame list of ZTF candidates (from Objects API)
    grb_date : GRB date in astropy.time format

    Returns
    -------
    list_pdf_ft : panda data frame list of the fast transient-selected candidates

    """
    
    list_pdf_ft = []
    
    for obj in list_pdf_objs:
        
        # add the mag dc and err_dc in the panda dat frame of the candidate
        # SHOULD BE INCLUDED BY DEFAULT IN FINK OBJECTS AT SOME POINT
        mag_dc, err_dc = np.transpose(
            [
                dc_mag(*args) for args in zip(
                    obj['i:fid'].astype(int).values,
                    obj['i:magpsf'].astype(float).values,
                    obj['i:sigmapsf'].astype(float).values,
                    obj['i:magnr'].astype(float).values,
                    obj['i:sigmagnr'].astype(float).values,
                    obj['i:magzpsci'].astype(float).values,
                    obj['i:isdiffpos'].values
                )
            ]
        )
        obj['i:magdc'] = mag_dc.tolist()
        obj['i:err_dc'] = err_dc.tolist()
        
        tag_ft = []

        for filt in np.unique(obj['i:fid']):
            mag_rate = []
            pdf_fid = get_pdf_measure(obj,filt)
            maskValid = pdf_fid['d:tag'] == 'valid'
            maskUpper= pdf_fid['d:tag'] == 'upperlim'
            pdf_fid = pdf_fid[maskValid | maskUpper].sort_values(by="i:jd")
    
            # if there is at least one measurement compute the pseudo mag rate 
            #between the first measured point and the previous U.L.
    
            if (pdf_fid['i:sigmapsf']>0).any():
                
                # get the index of the first measured data point
                index_meas = pdf_fid[pdf_fid['i:sigmapsf']>0]['i:magdc'].index[0]
                
                # get the index of the closest UL point prior the first measured one
                mask_ul_mag = pdf_fid['i:sigmapsf'].isnull()
                mask_ul_time = pdf_fid['i:jd']< pdf_fid['i:jd'][index_meas]
                index_ul = pdf_fid[mask_ul_mag & mask_ul_time]['i:magdc'].index[-1]
                
                # make the first pseudo mag-rate between the first data point 
                # and the last UL
                mag_meas = pdf_fid['i:magdc'][index_meas]
                time_meas = pdf_fid['i:jd'][index_meas]
                mag_ul = pdf_fid['i:diffmaglim'][index_ul]
                time_ul = pdf_fid['i:jd'][index_ul]
                
                # mag rate in mag/day (negative values = rising flux)
                mag_rate.append((mag_meas-mag_ul)/(time_meas-time_ul))
                
            if pdf_fid['v:rate(dg)'].notnull().any() and \
                (pdf_fid['v:rate(dg)']>0).any():
                mag_rate = mag_rate + pdf_fid[pdf_fid['v:rate(dg)'].\
                                              notnull()]['v:rate(dg)'].values.tolist()
    
            
            if (abs(np.array(mag_rate))>0.3).any():
                tag_ft.append(True)
            else:
                tag_ft.append(False)
        # BE CAREFUL, TO BE A FAST TRANSIENT WE REQUIRE TO HAVE A FT BEHAVIOR
        # (|dmag|>0.3mag/day) BOTH IN THE g AND R BANDS (CONSERVATIVE APPROACH)
        if np.array(tag_ft).all():
            list_pdf_ft.append(obj)
        
    return list_pdf_ft


def check_det_hist(pdf_cand,grb_date):
    """
    Find if a candidate have a detection history compatible with a GRB
    afterglow behavior

    Parameters
    ----------
    pdf_cand : pandas dataframe related to a ZTF object (from objects FINK API)

    Returns
    -------
    A True or False boolean

    """
    
    if pdf_cand['i:jd'] < grb_date.jd or pdf_cand['i:jd'] > grb_date.jd+90 or\
        pdf_cand['i:sigmapsf'] == 0:
        return True
    else:
        return False
