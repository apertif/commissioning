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
import aplpy

"""
Global variables
"""
mosaic_path = '/tank/apertif/mosaics/' #path to mosaic top level, assuming happili-05
internal_mosaic_path = 'mosaics/continuum/mosaic/' #internal path for continuum mosaics

def get_pybdsf_comp(taskid):
    """
    Basic function to take taskid and get pybdsf ouptput from continuum mosaic
    inputs:
        taskid: string of taskid
    outputs:
        data: table data from output csvfile
    """
    #get path
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

def get_beam_info(taskid):
    """
    Helper function to get beam sizes from header
    """
    filename = '{0}_mosaic.fits'.format(taskid)
    path_to_mosaic = os.path.join(mosaic_path,taskid,internal_mosaic_path,filename)
    hdr = getheader(path_to_mosaic)
    bmaj = hdr['BMAJ']*u.deg
    bmin = hdr['BMIN']*u.deg
    return bmaj.to(u.arcsec),bmin.to(u.arcsec)

def plot_maj_min(taskid):
    """
    Plot major vs. minor axes for sources
    """
    #first get data
    data = get_pybdsf_comp(taskid)
    
    #get beam sizes from header for comparison
    bmaj, bmin = get_beam_info(taskid)
    fig, (ax1, ax2) = plt.subplots(nrows=2, ncols=1, figsize = (8,12))
    
    ax1.scatter(data['Min']*u.deg.to(u.arcsec),data['Maj']*u.deg.to(u.arcsec),s=2)
    ax1.set_xlabel('Minor axis [arcsec]')
    ax1.set_ylabel('Major axis [arcsec]')
    
    ax1.set_xlim([0,60])
    ax1.set_ylim([0,60])
    
    ax1.plot([0,60],[bmaj.value,bmaj.value],'k:')
    ax1.plot([bmin.value,bmin.value],[0,60],'k:')
    ax1.plot([0,60],[0,60],'k')
    
    ax2.scatter(data['Peak_flux']*u.Jy.to(u.mJy),data['Min']/data['Maj'],s=2)
    ax2.set_ylabel('Minor / Major')
    ax2.set_xlabel('Peak flux (mJy/bm)')
    ax2.set_xlim([0,100])
    
    return fig
    
def plot_pa(taskid):
    """
    Plot PA for sources
    """
    #first get data
    data = get_pybdsf_comp(taskid)
    
    #get beam sizes from header for comparison
    fig, (ax1) = plt.subplots(nrows=1, ncols=1, figsize = (8,6))
    
    ax1.scatter(data['Min']/data['Maj'],data['PA'],s=2)
    ax1.set_xlabel('Minor / Major')
    ax1.set_ylabel('PA [deg]')
   
    return fig

def plot_peak_size(taskid,image_plane=False):
    """
    Plot peak flux against major and minor axes
    """
    #first get data
    data = get_pybdsf_comp(taskid)
    
    
    #set up figure
    fig, (ax1,ax2) = plt.subplots(nrows=2, ncols=1, figsize = (8,12), 
                                  sharey = True)
    
    
    if image_plane == True:
        ax1.scatter(data['Maj_img_plane']*u.deg.to(u.arcsec),data['Peak_flux']*u.Jy.to(u.mJy),s=3)
        ax2.scatter(data['Min_img_plane']*u.deg.to(u.arcsec),data['Peak_flux']*u.Jy.to(u.mJy),s=3)
        ax1.set_xlabel('Major axis, image plane [arcsec]')
        ax2.set_xlabel('Minor axis, image plane [arcsec]')
    else:
        ax1.scatter(data['Maj']*u.deg.to(u.arcsec),data['Peak_flux']*u.Jy.to(u.mJy),s=3)
        ax2.scatter(data['Min']*u.deg.to(u.arcsec),data['Peak_flux']*u.Jy.to(u.mJy),s=3)
        ax1.set_xlabel('Major axis [arcsec]')
        ax2.set_xlabel('Minor axis [arcsec]')  
    
    ax1.set_ylabel(' Peak (mJy/bm) ')
    
    ax2.set_ylabel(' Peak (mJy/bm) ')    
    
    bmaj, bmin = get_beam_info(taskid)
   
    ax1.set_ylim([0,500])
    ax2.set_ylim([0,500])
    
    ax1.plot([bmaj.value,bmaj.value],[0,1000],'k:')
    ax2.plot([bmin.value,bmin.value],[0,1000],'k:')
    
    return fig

def plot_total_size(taskid):
    """
    Plot total flux against major and minor axes
    """
    #first get data
    data = get_pybdsf_comp(taskid)
    
    
    #set up figure
    fig, (ax1,ax2) = plt.subplots(nrows=2, ncols=1, figsize = (8,12), 
                                  sharey = True)
    
    ax1.scatter(data['Maj']*u.deg.to(u.arcsec),data['Total_flux']*u.Jy.to(u.mJy),s=3)
    ax1.set_xlabel('Major axis [arcsec]')
    ax1.set_ylabel(' Total (mJy/bm) ')
    
    ax2.scatter(data['Min']*u.deg.to(u.arcsec),data['Total_flux']*u.Jy.to(u.mJy),s=3)
    ax2.set_xlabel('Minor axis [arcsec]')
    ax2.set_ylabel(' Total (mJy/bm) ')
    
    
    bmaj, bmin = get_beam_info(taskid)
    
    ax1.set_ylim([0,500])
    ax2.set_ylim([0,500])
   
    ax1.plot([bmaj.value,bmaj.value],[0,1000],'k:')
    ax2.plot([bmin.value,bmin.value],[0,1000],'k:')
    
    return fig


def plot_cutouts_minor(taskid, minormin, minormax):
    """
    Plot cutouts of sources with minor axis limits as given
    inputs:
        taskid: string
        minormin: smallest minor axis size in arcsec (units)
        minormax: largest minor axis in arcsec (units)
    """
    filename = '{0}_mosaic.fits'.format(taskid)
    path_to_mosaic = os.path.join(mosaic_path,taskid,internal_mosaic_path,filename)
    
    data = get_pybdsf_comp(taskid)
    
    #find indices of sources
    minoravg = (minormin + minormax) / 2. 
    minorhalf = (minormax - minormin) / 2.
    ind = np.where(np.abs(data['Min']-minoravg.to(u.deg).value) <= 
                   minorhalf.to(u.deg).value)[0]

    nplots = len(ind)
    nx = 4.
    ny = int(np.ceil(nplots / nx))
        
    fig = plt.figure(figsize=(3*nx, 3*ny))
    
    if nplots > 20:
        print(('Requesting {0} plots! This is too many. '
               'Please retry minor axis limits').format(nplots))
    else:  
        print('Plotting {0} sources'.format(nplots))
        i=1
 
        for source in ind:
        #iterate over index array
            try:
                print(i)
                f1 = aplpy.FITSFigure(path_to_mosaic,figure=fig,subplot=(ny,nx,i))
                f1.show_grayscale()
                #recenter w/ 2 arcmin radius
                f1.recenter(data['RA'][source],data['DEC'][source],radius=2/60.)
                #plot fitted shape
                f1.show_ellipses(data['RA'][source],data['DEC'][source],
                                 data['Min'][source], #15/3600., #width
                                 data['Maj'][source], #30/3600., #height
                                 angle=data['PA'][source], #angle=0, #angle
                                 edgecolor='red',facecolor='None')
                i=i+1
            except:
                #if a plot fails, keep going
                i=i+1
        plt.show();
                
    return fig
    