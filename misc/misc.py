#python modules for miscellaneous commissioning tests

#load necessary packages
import os
import numpy as np
import sys
sys.path.append('/home/adams/altadata')
sys.path.append('/home/adams/commissioning/crosscal')
import crosscal as cc
import casacore.tables as pt
from astropy.coordinates import SkyCoord
import matplotlib.pyplot as plt

"""
First set of miscellaneous tests is to look at the central positions of different beams
when they are specified to have the target in them. This appears to not be 100% correct
which has implications for coordinates and phase tracking of outlying beams.
This will build on some of the functionality in crosscal for specifying a set of switching
scans and retrieving the data. It will add retrieving all data for the beam00 pointing for
reference positions. It will also add the retrieval of coordinates, calculation of offsets,
and visualization of the offsets.
"""

def get_data_beam_post(startscan,endscan,refscan,obsrecordfile):
    #Get data for a startscan,endscan, and a scan specified as reference to get relative positions
    #For the switching observations, use crosscal infrastructure:
    scans = cc.ScanSpecification()
    scans.setstartscan(startscan)
    scans.setendscan(endscan)
    cc.copy_scans(scans,obsrecordfile,'/data/adams/apertif/beampos',run=True)
    
    #for the reference scan, want to retrieve all beams
    #This is so I can have coordinates/relative positions of beams for visualization
    #Since I changed directory above, I should still be there (hope!)
    altadata_string_command = "python /home/adams/altadata/getdata_alta.py {0}-{0} 00-36".format(refscan)
    os.system(altadata_string_command)
    
def get_coords_and_offsets(startscan,endscan,refscan):
    #want to get my coordinates and offsets for full range of scans
    #first, get scan & beam list to help make life easier:
    scans = cc.ScanSpecification()
    scans.setstartscan(startscan)
    scans.setendscan(endscan)
    mode,scan_list,beam_list = cc.get_scan_list(scans,obsrecordfile)
    
    #then set arrays of the right length to store coordinates
    coord_list = np.empty(len(scan_list),dtype=object)
    offset_list = np.empty(len(scan_list),dtype=object)
    #get reference coordinate
    refMS = 'WSRTA{0}_B000.MS'.format(refscan)
    t_field = pt.table(msfile+"::FIELD", readonly=True)
    phasedir=t_field[0]["PHASE_DIR"]
    ref_coord = SkyCoord(phasedir[0],phasedir[1],unit='rad')
    #now go through each observation to get coordinates
    for i,scan,beam in enumerate(zip(scan_list,beam_list)):
        #format the msfile:
        msfile = 'WSRTA{0}_B{1:0>3}.MS'.format(scan,beam)
        #read the FIELD table
        t_field = pt.table(msfile+"::FIELD", readonly=True)
        #take PHASE_DIR as center coordinate of beam
        phasedir=t_field[0]["PHASE_DIR"]
        #this is [ra,dec] in radians
        #put into a SkyCoord object - easy to handle
        c=SkyCoord(phasedir[0],phasedir[1],unit='rad')
        coord_list[i] = c
        offset_list[i] = c.separation(ref_coord).value
    return coord_list,offset_list,beam_list

def print_offsets(startscan,endscan,refscan):
    #produce a nicely formatted overview of offsets
    coord_list,offset_list,beam_list = get_coords_and_offsets(startscan,endscan,refscan)
    print 'Beam  Offset in degrees'
    for beam,offset in zip(beam_list,offset_list):
        print 'B{0:0>3}  {1}'.format(beam,offset)

"""
Here I may want to define a function for doing the nice visualization
But probably best/easiest to start off playing with this in a notebook
def visualize_offsets(startscan,endscan,refscan):
    #do nice visualization of the offset as function of position
    #to do this, need to get reference positions for every beam along the way

    #start with basic information
    coord_list,offset_list,beam_list = get_coords_and_offsets(startscan,endscan,refscan)
    
    #initialize visualization plot
    #get reference direction    
    refMS = 'WSRTA{0}_B000.MS'.format(refscan)
    t_field = pt.table(msfile+"::FIELD", readonly=True)
    phasedir=t_field[0]["PHASE_DIR"]
    ref_coord = SkyCoord(phasedir[0],phasedir[1],unit='rad')
    fig = aplpy.FITSFigure(rcimage,figure=fig,aspect='equal'
    
    #iterate through beam
    #plot each, and get reference position as part of this
    for beam,offset,coord in zip(beam_list,offset_list,coord_list):
        
    return fig
 
 """
    
        