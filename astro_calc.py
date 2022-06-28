#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  4 11:04:48 2021

@author: Damien Turpin : damien.turpin@cea.fr
"""

from astropy.coordinates import SkyCoord
import numpy as np
from math import pi, sqrt
from scipy.stats import poisson
from scipy import special
from scipy.optimize import minimize

def ang_distance(coord_grb,coord_ztf):
    c1 = SkyCoord(coord_grb[0], coord_grb[1], unit='degree', frame='icrs')
    if type(coord_ztf[0]) == list:
        sep = []
        for coord in coord_ztf:
            c2 = SkyCoord(coord[0], coord[1], unit='degree', frame='icrs')
            sep.append(c1.separation(c2))
    else:
        c2 = SkyCoord(coord_ztf[0], coord_ztf[1], unit='degree', frame='icrs')
        sep = c1.separation(c2)
    
    
    return sep

def p_ser_grb(error_radius,size_time_window, r_grb):
    """
    Created on Mon Oct  4 10:34:09 2021
    
    @author: Damien Turpin : damien.turpin@cea.fr
    
    function that gives the chance probability of having a positive spatial and 
    temporal match between a GRB and a ZTF transient candidate
    Inputs:
        error radius of the GRB localization region in degree: error_radius
        size of the searching time window in year: size_time_window
        GRB detection rate for a given satellite in events/year : R_GRB
    Outputs:
        Serendipituous probabilities for a GRB/ZTF candidate association
    """

    # omega = 2*pi*(1-cos(radians(error_radius))) # solid angle in steradians
    grb_loc_area = pi*(error_radius)**2 # in square degrees
    allsky_area = 4*pi*(180/pi)**2 # in square degrees
    ztf_coverage_rate = 3750 # sky coverage rate of ZTF in square degrees per hour
    limit_survey_time = 4 #duration (in hour) during which ZTF will cover individual parts of the sky in a night
    
    # short and long GRB detection rate
    r_sgrb = r_grb/3
    r_lgrb = r_grb-r_sgrb
    
    # Poisson probability of detecting a GRB during a searching time window
    
    p_grb_detect_ser = 1-poisson.cdf(1,r_grb*size_time_window)
    p_lgrb_detect_ser = 1-poisson.cdf(1,r_lgrb*size_time_window)
    p_sgrb_detect_ser = 1-poisson.cdf(1,r_sgrb*size_time_window)
    
    # we limit the fraction of the sky ZTF is able to cover to 4 hours of continuous survey
    # we consider that every day (during several days only) ZTF will cover the same part of
    # the sky with individual shots (so no revisit) during 4 hours
        
    if size_time_window*365.25*24 <= limit_survey_time:
        ztf_sky_frac_area = (ztf_coverage_rate*size_time_window*365.25*24)
    else:
        ztf_sky_frac_area = ztf_coverage_rate*limit_survey_time
    
    # probability of finding a GRB within the region area paved by ZTF during a given amount of time
    p_grb_in_ztf_survey = (ztf_sky_frac_area/allsky_area) * p_grb_detect_ser
    p_lgrb_in_ztf_survey = (ztf_sky_frac_area/allsky_area) * p_lgrb_detect_ser
    p_sgrb_in_ztf_survey = (ztf_sky_frac_area/allsky_area) * p_sgrb_detect_ser
    
    # probability of finding a ZTF transient candidate inside the GRB error box 
    # knowing the GRB is in the region area paved by ZTF during a given amount of time

    p_ser_grb = p_grb_in_ztf_survey*(grb_loc_area/ztf_sky_frac_area)
    
    p_ser_lgrb = p_lgrb_in_ztf_survey*(grb_loc_area/ztf_sky_frac_area)
    
    p_ser_sgrb = p_sgrb_in_ztf_survey*(grb_loc_area/ztf_sky_frac_area)
    
    p_sers = [p_ser_grb, p_ser_lgrb, p_ser_sgrb]
    
    return p_sers

def sig_est(prob):
    fun = lambda x: abs(prob-special.erf(x/sqrt(2)))
    res = minimize(fun, [0] , method='Nelder-Mead')
    
    return res.x

def mag2fluxcal_snana(magpsf: float, sigmapsf: float):
    """ Conversion from magnitude to Fluxcal from SNANA manual
    Parameters
    ----------
    magpsf: float
        PSF-fit magnitude from ZTF
    sigmapsf: float
    Returns
    ----------
    fluxcal: float
        Flux cal as used by SNANA
    fluxcal_err: float
        Absolute error on fluxcal (the derivative has a minus sign)
    """
    if magpsf is None:
        return None, None
    fluxcal = 10 ** (-0.4 * magpsf) * 10 ** (11)
    fluxcal_err = 9.21034 * 10 ** 10 * np.exp(-0.921034 * magpsf) * sigmapsf

    return fluxcal, fluxcal_err

def apparent_flux(fid, magpsf, sigmapsf, magnr, sigmagnr, magzpsci, isdiffpos):
    """ Compute apparent flux from difference magnitude supplied by ZTF
    This was heavily influenced by the computation provided by Lasair:
    https://github.com/lsst-uk/lasair/blob/master/src/alert_stream_ztf/common/mag.py
    Paramters
    ---------
    fid
        filter, 1 for green and 2 for red
    magpsf,sigmapsf; floats
        magnitude from PSF-fit photometry, and 1-sigma error
    magnr,sigmagnr: floats
        magnitude of nearest source in reference image PSF-catalog
        within 30 arcsec and 1-sigma error
    magzpsci: float
        Magnitude zero point for photometry estimates
    isdiffpos: str
        t or 1 => candidate is from positive (sci minus ref) subtraction;
        f or 0 => candidate is from negative (ref minus sci) subtraction
    Returns
    --------
    dc_flux: float
        Apparent magnitude
    dc_sigflux: float
        Error on apparent magnitude
    """
    if magpsf is None:
        return None, None
    # zero points. Looks like they are fixed.
    ref_zps = {1: 26.325, 2: 26.275, 3: 25.660}
    magzpref = ref_zps[fid]

    # reference flux and its error
    magdiff = magzpref - magnr
    if magdiff > 12.0:
        magdiff = 12.0
    ref_flux = 10**(0.4 * magdiff)
    ref_sigflux = (sigmagnr / 1.0857) * ref_flux

    # difference flux and its error
    if magzpsci == 0.0:
        magzpsci = magzpref
    magdiff = magzpsci - magpsf
    if magdiff > 12.0:
        magdiff = 12.0
    difference_flux = 10**(0.4 * magdiff)
    difference_sigflux = (sigmapsf / 1.0857) * difference_flux

    # add or subract difference flux based on isdiffpos
    if isdiffpos == 't':
        dc_flux = ref_flux + difference_flux
    else:
        dc_flux = ref_flux - difference_flux

    # assumes errors are independent. Maybe too conservative.
    dc_sigflux = np.sqrt(difference_sigflux**2 + ref_sigflux**2)

    return dc_flux, dc_sigflux

def dc_mag(fid, magpsf, sigmapsf, magnr, sigmagnr, magzpsci, isdiffpos):
    """ Compute apparent magnitude from difference magnitude supplied by ZTF
    Parameters
    Stolen from Lasair.
    ----------
    fid
        filter, 1 for green and 2 for red
    magpsf,sigmapsf
        magnitude from PSF-fit photometry, and 1-sigma error
    magnr,sigmagnr
        magnitude of nearest source in reference image PSF-catalog
        within 30 arcsec and 1-sigma error
    magzpsci
        Magnitude zero point for photometry estimates
    isdiffpos
        t or 1 => candidate is from positive (sci minus ref) subtraction;
        f or 0 => candidate is from negative (ref minus sci) subtraction
    """
    # zero points. Looks like they are fixed.
    ref_zps = {1: 26.325, 2: 26.275, 3: 25.660}
    magzpref = ref_zps[fid]

    # difference flux and its error
    if magzpsci is None:
        magzpsci = magzpref

    dc_flux, dc_sigflux = apparent_flux(
        fid, magpsf, sigmapsf, magnr, sigmagnr, magzpsci, isdiffpos
    )

    # apparent mag and its error from fluxes
    if (dc_flux == dc_flux) and dc_flux > 0.0:
        dc_mag = magzpsci - 2.5 * np.log10(dc_flux)
        dc_sigmag = dc_sigflux / dc_flux * 1.0857
    else:
        dc_mag = magzpsci
        dc_sigmag = sigmapsf

    return dc_mag, dc_sigmag