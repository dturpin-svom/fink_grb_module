#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 24 11:05:42 2022

@author: dt261638
"""

import matplotlib.pyplot as plt
import numpy as np

def plot_cand_lc(pdf, grb_date,size_time_window):
    """

    Parameters
    ----------
    pdf : object panda dataframe related to one ZTF object 
        DESCRIPTION.

    Returns
    -------
    None.

    """
    # Plot the light curves of the fast rising/fading transients
    colordic = {1: 'C0', 2: 'C1'}
    filter_label = {1: 'g', 2: 'R'}
    
    
    # filter out the data point before the GRB  trigger time
    mask_grbtime = pdf['i:jd']>grb_date.jd-30
    # setup the time unit for the light curve
    unit_time_factor = 1 # here in hour (86400 in sec, 1440 in min, 24 in hour)
    if unit_time_factor == 1:
        unit_time = 'days'
    elif unit_time_factor == 24:
        unit_time = 'hours'
    elif unit_time_factor == 1440:
        unit_time = 'min'
    elif unit_time_factor == 86400:
        unit_time = 'sec'
    for filt in np.unique(pdf['i:fid']):
        fig = plt.figure(figsize=(15, 6))
        maskFilt = pdf['i:fid'] == filt
            
        # The column `d:tag` is used to check data type
        maskValid = pdf['d:tag'] == 'valid'
        plt.errorbar(
            pdf[maskValid & maskFilt & mask_grbtime]['i:jd'].\
            apply(lambda x: (x - grb_date.jd)*unit_time_factor),
            pdf[maskValid & maskFilt & mask_grbtime]['i:magdc'],
            pdf[maskValid & maskFilt & mask_grbtime]['i:err_dc'],
            ls = '--', marker='o', color=colordic[filt], label=['valid']
        )
    
        maskUpper = pdf['d:tag'] == 'upperlim'
        plt.plot(
            pdf[maskUpper & maskFilt & mask_grbtime]['i:jd'].\
            apply(lambda x: (x - grb_date.jd)*unit_time_factor),
            pdf[maskUpper & maskFilt & mask_grbtime]['i:diffmaglim'],
            ls='', marker='^', color=colordic[filt], markerfacecolor='none', label=['upperlim']
        )
    
        maskBadquality = pdf['d:tag'] == 'badquality'
        plt.errorbar(
            pdf[maskBadquality & maskFilt & mask_grbtime]['i:jd'].\
            apply(lambda x: (x - grb_date.jd)*unit_time_factor),
            pdf[maskBadquality & maskFilt & mask_grbtime]['i:magdc'],
            pdf[maskBadquality & maskFilt & mask_grbtime]['i:err_dc'],
            ls='', marker='v', color=colordic[filt], label=['badquality']
        )
        if pdf[maskValid & maskFilt]['i:magpsf'].size > 0:
            plt.plot([0,0],
                     [min(pdf[maskValid & maskFilt & mask_grbtime]['i:magdc']),20.5],
                     ls='--',lw=3,c='red',label='GRB trigger time')
            plt.plot([size_time_window*unit_time_factor,size_time_window*unit_time_factor],
                     [min(pdf[maskValid & maskFilt & mask_grbtime]['i:magdc']),20.5],
                     ls='--',lw=3,c='black',label='End of searching time window')
    
        plt.gca().invert_yaxis()
    #         plt.xscale('log')
        plt.xlabel('Time since the GRB trigger time ['+unit_time+']')
        plt.ylabel('Magnitude')
        plt.title("ZTF light curve of "+np.unique(pdf['i:objectId'].values)+" in "+\
                  filter_label[filt]+"-band")
        plt.legend(loc='best')
        plt.show()