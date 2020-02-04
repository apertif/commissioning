#Helper functions for looking at flux issues

from __future__ import print_function

__author__ = "E.A.K. Adams"

"""
Helper functions for making figures
for SVC data release paper.
First set of functions is for plotting
data for discussion of flagging.
"""

import matplotlib.pyplot as plt
import os
import numpy as np
import astropy.units as u
from astropy.io.fits import getheader
from astropy.io import ascii

def get_pybdsf_comp(taskid):
    """
    Basic function to take taskid and get pybdsf ouptput from continuum mosaic
    inputs:
        taskid: string of taskid
    outputs:
        data: table data from output csvfile
    """
    #get path
    mosaic_path = '/tank/apertif/mosaics/'
    internal_mosaic_path = 'mosaics/continuum/mosaic/'
    filename = '{0}_mosaic_pybdsf_comp.csv'.format(taskid)
    path_to_comp = os.path.join(mosaic_path,taskid,internal_mosaic_path,filename)
    #load file
    data = ascii.read(path_to_comp,header_start=4,format='csv',data_start=5)
    return data

def plot_int_peak_norm_size(taskid):
    """
    Plot integrated over peak flux vs. major and minor axis (two plots)
    inputs:
        taskid - string w/ taskid of continuum mosaic
    """
    #first get data
    data = get_pybdsf_comp(taskid)
    
    rel_flux_diff = (data['Total_flux'] - data['Peak_flux']) / data['Peak_flux']
    
    #set up figure
    fig, (ax1,ax2) = plt.subplots(nrows=2, ncols=1, figsize = (8,12), 
                                  sharey = True)
    
    ax1.scatter(data['Maj']*u.deg.to(u.arcsec),rel_flux_diff)
    ax1.set_xlabel('Major axis [arcsec]')
    ax1.set_ylabel(' (Total - Peak) / Peak ')
    
    ax2.scatter(data['Min']*u.deg.to(u.arcsec),rel_flux_diff)
    ax2.set_xlabel('Minor axis [arcsec]')
    ax2.set_ylabel(' (Total - Peak) / Peak ')
    
    ax1.set_xlim([0,60])
    ax2.set_xlim([0,60])
    
    ax1.plot([0,60],[0,0],'k:')
    ax2.plot([0,60],[0,0],'k:')
    
    return fig

def plot_int_peak_size(taskid):
    """
    Plot integrated over peak flux vs. major and minor axis (two plots)
    inputs:
        taskid - string w/ taskid of continuum mosaic
    """
    #first get data
    data = get_pybdsf_comp(taskid)
    
    flux_ratio = data['Total_flux'] / data['Peak_flux']
    
    #set up figure
    fig, (ax1,ax2) = plt.subplots(nrows=2, ncols=1, figsize = (8,12), 
                                  sharey = True)
    
    ax1.scatter(data['Maj']*u.deg.to(u.arcsec),flux_ratio)
    ax1.set_xlabel('Major axis [arcsec]')
    ax1.set_ylabel(' Total / Peak ')
    
    ax2.scatter(data['Min']*u.deg.to(u.arcsec),flux_ratio)
    ax2.set_xlabel('Minor axis [arcsec]')
    ax2.set_ylabel(' Total / Peak ')
    
    ax1.set_xlim([0,60])
    ax2.set_xlim([0,60])
    
    ax1.plot([0,60],[1,1],'k:')
    ax2.plot([0,60],[1,1],'k:')
    
    return fig

def plot_maj_min(taskid):
    """
    Plot major vs. minor axes for sources
    """
    #first get data
    data = get_pybdsf_comp(taskid)
    
    #get beam sizes from header for comparison
    bmaj, bmin = get_beam_info('191209026')
    fig, (ax1) = plt.subplots(nrows=1, ncols=1, figsize = (8,6))
    
    ax1.scatter(data['Min']*u.deg.to(u.arcsec),data['Maj']*u.deg.to(u.arcsec))
    ax1.set_xlabel('Minor axis [arcsec]')
    ax1.set_ylabel('Major axis [arcsec]')
    
    ax1.set_xlim([0,60])
    ax1.set_ylim([0,60])
    
    ax1.plot([0,60],[bmaj.value,bmaj.value],'k:')
    ax1.plot([bmin.value,bmin.value],[0,60],'k:')
    ax1.plot([0,60],[0,60],'k')
    
    return fig
    