#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  4 14:04:48 2021

@author: Damien Turpin : damien.turpin@cea.fr
"""

import pandas as pd
from astropy.table import Table
import numpy as np
from astropy.time import Time


def convert_time(time,format_time):
    
    return Time(time,format = format_time)

def cat_select(mission_name,cat_path):
    
    # fermi_cat = Table.read(cat_file)
    if mission_name == 'Fermi':
        delimiter = '|'
        cat_table = Table.read(cat_path)
        pdf = pd.DataFrame(np.array(cat_table))
        pdf['trigger_time'] = pdf['trigger_time'].\
            apply(lambda x: convert_time(x,'mjd'))
    elif mission_name == 'Swift':
        delimiter = ','
        pdf = pd.read_csv(cat_path, sep = delimiter)
        pdf['trigger_time'] = pdf['trigger_time'].\
            apply(lambda x: convert_time(x,'iso'))
        pdf['error_radius'] = pdf['error_radius'].\
            apply(lambda x: x/60)
    return pdf



