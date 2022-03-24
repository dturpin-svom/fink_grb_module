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
        'startdate_conesearch': grb_date.iso, # optional
        'window_days_conesearch': size_time_window, # search for ZTF candidates up to 7 days post GRB
      }
    )
    
    # Format output in a DataFrame
    return pd.read_json(r.content)


def get_fink_objects(ztf_objectIds):
    
    # get data for many objects
    r = requests.post(
      'https://fink-portal.org/api/v1/objects',
      json={
        'objectId': ','.join(ztf_objectIds),
        'output-format': 'json',
        'withupperlim': 'True'
      }
    )
    
    # Format output in a DataFrame
    pdf = pd.read_json(r.content)
        
    list_pdf_objs = []
    dict_pdf_objs = {}
    for obj in ztf_objectIds:
        dict_pdf_objs[obj] = pdf[pdf['i:objectId']==obj]
        list_pdf_objs.append(pdf[pdf['i:objectId']==obj])
    
    return list_pdf_objs, dict_pdf_objs