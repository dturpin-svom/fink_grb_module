#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  4 11:14:48 2021

@author: Damien Turpin : damien.turpin@cea.fr

A set of routines to query the Fink science portal

"""

import requests

import pandas as pd
import numpy as np

def ft_finder(pdf,grb_date):
    """

    Search for fast-fading transients in the positive GRB/ZTF crossmatched list
    ----------
    pdf : panda data frame of a list of ZTF candidates
    grb_date : GRB date in astropy.time format

    Returns
    -------
    pdf_ft : panda data frame of the fast transient-selected candidates

    """

    for obj in pdf['i:objectId'].values:


        # get data for the candidate
        r = requests.post(
          'https://fink-portal.org/api/v1/objects',
          json={
            'objectId': obj,
            'output-format': 'json',
            'withupperlim': 'True'
          }
        )


        # Format output in a DataFrame
        pdf_cand = pd.read_json(r.content)
        
        # we do not keep transient candidates with detections previous to the GRB trigger time
        nodethist_filter = pdf_cand.loc[pdf_cand['i:jd'] < grb_date.jd,
                                        'i:sigmapsf'].isnull().values.all()
        # we do not keep transient candidates with detections 3 month post GRB trigger time
        nodet_postgrb_filter = pdf_cand.loc[pdf_cand['i:jd'] > grb_date.jd+90,
                                            'i:sigmapsf'].isnull().values.all()
        # ensure we do not keep the transients with a detection history prior to the GRB trigger time
        if nodethist_filter == True and nodet_postgrb_filter == True:

            for filt in np.unique(pdf_cand['i:fid']):
                maskFilt = pdf_cand['i:fid'] == filt

                if filt == 1:
                    maskValid = pdf_cand['d:tag'] == 'valid'
                    maskUpper= pdf_cand['d:tag'] == 'upperlim'
                    pdf_cand_gband = pdf_cand[maskFilt & \
                                              (maskValid | maskUpper)].sort_values(by="i:jd")
                elif filt == 2:
                    maskValid = pdf_cand['d:tag'] == 'valid'
                    maskUpper= pdf_cand['d:tag'] == 'upperlim'
                    pdf_cand_rband = pdf_cand[maskFilt & \
                                              (maskValid | maskUpper)].sort_values(by="i:jd")

            # compute the g-band rising/fading rate as function of time (récupérer fonction J.Peloton pour calcul apparent mag ou DC mag)
            mag_rate_gband = []
            if len(pdf_cand_gband[pdf_cand_gband['i:sigmapsf']>0])>1:
                mag_diff = pdf_cand_gband['i:magpsf']\
                    [pdf_cand_gband['i:sigmapsf'].isnull()==False].values[0]-\
                        pdf_cand_gband['i:diffmaglim'].values[-1]

                time_diff = abs(pdf_cand_gband['i:jd']\
                                [pdf_cand_gband['i:sigmapsf'].isnull()].values[-1]-\
                                    pdf_cand_gband['i:jd']\
                                        [pdf_cand_gband['i:sigmapsf'].isnull()==False].values[0])

                mag_rate_gband .append(mag_diff/time_diff)
                k = 1
                while k<len(pdf_cand_gband['v:rate(dg)'][pdf_cand_gband['i:sigmapsf']>0]) :
                    mag_rate = pdf_cand_gband['v:rate(dg)'][pdf_cand_gband['i:sigmapsf']>0].values[k]
                    mag_rate_gband.append(mag_rate)
                    k = k+1

            # compute the R-band rising/fading rate as function of time   
            mag_rate_rband = []
            if len(pdf_cand_rband[pdf_cand_rband['i:sigmapsf']>0])>1:
                mag_diff = pdf_cand_rband['i:magpsf']\
                    [pdf_cand_rband['i:sigmapsf'].isnull()==False].values[0]-\
                        pdf_cand_rband['i:diffmaglim'].values[-1]

                time_diff = abs(pdf_cand_rband['i:jd']\
                                [pdf_cand_rband['i:sigmapsf'].isnull()].values[-1]-\
                                    pdf_cand_rband['i:jd']\
                                        [pdf_cand_rband['i:sigmapsf'].isnull()==False].values[0])

                mag_rate_rband .append(mag_diff/time_diff)
                k = 1
                while k<len(pdf_cand_rband['v:rate(dr)'][pdf_cand_rband['i:sigmapsf']>0]) :
                    mag_rate = pdf_cand_rband['v:rate(dr)'][pdf_cand_rband['i:sigmapsf']>0].values[k]
                    mag_rate_rband.append(mag_rate)
                    k = k+1
    #         elif len(pdf_cand_rband[pdf_cand_rband['i:sigmapsf']>0]) == 1:
    #             mag_diff = pdf_cand_rband['i:magpsf'][pdf_cand_rband['i:sigmapsf'].isnull()==False].values[0]-\
    #             pdf_cand_rband['i:diffmaglim'].values[-1]


    #             time_diff = abs(pdf_cand_rband['i:jd'][pdf_cand_rband['i:sigmapsf'].isnull()].values[-1]-\
    #             pdf_cand_rband['i:jd'][pdf_cand_rband['i:sigmapsf'].isnull()==False].values[0])

    #             mag_rate_rband.append(mag_diff/time_diff)


            # tag the fast rising or fading transients

            if len([i for i in mag_rate_gband if abs(i) >= 0.3])>0:
                gband_fast_tag = 1
            else:
                gband_fast_tag = 0

            if len([i for i in mag_rate_rband if abs(i) >= 0.3])>0:
                rband_fast_tag = 1
            else:
                rband_fast_tag = 0

            # Select only the candidates that have at least one fast evolving tag
            if rband_fast_tag == 1 or gband_fast_tag ==1:
                pdf.loc[pdf['i:objectId'] == obj,'v:ft'] = True

    return pdf.loc[pdf['v:ft'] == True]


def check_det_hist(objectId,grb_date):
    """
    Find if a candidate have a detection history compatible with a GRB
    afterglow behavior

    Parameters
    ----------
    objectId : ZTF objectId in string
    grb_date : GRB date in astropy.time format

    Returns
    -------
    A True or False boolean

    """
    
    # get data for the candidate
    r = requests.post(
      'https://fink-portal.org/api/v1/objects',
      json={
        'objectId': objectId,
        'output-format': 'json',
        'withupperlim': 'True'
      }
    )
    # Format output in a DataFrame
    pdf_cand = pd.read_json(r.content)
    
    # we do not keep transient candidates with detections previous to the GRB trigger time
    nodethist_filter = pdf_cand.loc[pdf_cand['i:jd'] < grb_date.jd,
                                    'i:sigmapsf'].isnull().values.all()
    # we do not keep transient candidates with detections 3 month post GRB trigger time
    nodet_postgrb_filter = pdf_cand.loc[pdf_cand['i:jd'] > grb_date.jd+90,
                                        'i:sigmapsf'].isnull().values.all()
    
    if nodethist_filter and nodet_postgrb_filter:
        return True
    else:
        return False
    
def explorer_crossmatch(ra,dec,cone_radius,grb_date,size_time_window):
    
    """
    Time and spatial crossmatch between the Fermi GRB catalog 
    and the ZTF transient candidate
    
    """
    r = requests.post(
      'https://fink-portal.org/api/v1/explorer',
      json={
        'ra': ra,
        'dec': dec,
        'radius': cone_radius,
        'startdate_conesearch': grb_date.isot, # optional
        'window_days_conesearch': size_time_window, # search for ZTF candidates up to 7 days post GRB
      }
    )
    
    # Format output in a DataFrame
    return pd.read_json(r.content)


    
    