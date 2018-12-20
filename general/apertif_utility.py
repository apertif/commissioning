#General utility module for Apertif
#Use for functions and such that want for multiple purposes

#load necessary packages
import os
import numpy as np
from astropy.io import ascii
import apercal
import casacore.tables as pt

def format_phasedir(phasedir):
    #takes phase_dir from measurement set and 
    #format it into HH:MM:SS DD:MM:SS
    ra_deg = phasedir[0,0] * 180./np.pi
    if ra_deg <0:
        ra_deg = ra_deg+360.
    dec_deg = phasedir[0,1]*180./np.pi
    rahr =int( np.floor(ra_deg/15.))
    ram = int(np.floor((ra_deg/15.-rahr)*60.))
    ras = ((ra_deg/15.-rahr)*60. - ram)*60.
    decd = int(np.floor(dec_deg))
    decm = int(np.floor( (dec_deg - decd)*60. ))
    decs =int( ( (dec_deg - decd)*60.  - decm )*60.)
    string = '{0:d}:{1:d}:{2:.1f} {3:d}:{4:d}:{5:d}'.format(rahr,ram,ras,decd,decm,decs)
    return string

def flip_ra(msfile):
    #Take a measurement set and flip RA about central pointing (Ref)
    #this takes full MS path
    #Iterating over beam would happen on a different level
    t_field = pt.table(msfile+"::FIELD", readonly=False)  #open MS for writing
    phasedir=t_field[0]["PHASE_DIR"]  #get the phasedir
    delaydir = t_field[0]["DELAY_DIR"] #get the delaydir, what I actually want to use for updating phase_dir
    refdir = t_field[0]["REFERENCE_DIR"] #get the reference dir, to reflect phasedir around
    newphasedir = np.copy(phasedir)
    newphasedir[0,0] = delaydir[0,0]+2*(refdir[0,0]-delaydir[0,0])  #reflect the delaydir around refdir to get new phasedir coord
    #This means that can run this code on a MS multiple times and it will work
    string_delay = format_phasedir(delaydir)
    string_phase = format_phasedir(newphasedir)
    print 'The delay direction is {0} and the new phase direction is {1}'.format(string_delay,string_phase)
    print ''
    #print to check against delay direction
    
    #do the actual update
    t_field.putcell("PHASE_DIR", 0, newphasedir)  #update the phase direction
    t_field.flush() #make sure changes are written to MS
    
    #recalculate uv coordinates
    pt.taql('update {0} set UVW = mscal.newuvw()'.format(msfile)) 
    
    
def unflip_ra(msfile):
    #In case accidentally flip something that shouldn't be, want to flip it back
    #This means revert to delay_dir
    #Note that if we ever change things so that delay_dir and phase_dir are not the same
    #This will (desperately) need to be updated
    t_field = pt.table(msfile+"::FIELD", readonly=False)  #open MS for writing
    phasedir=t_field[0]["PHASE_DIR"]  #get the phasedir
    delaydir = t_field[0]["DELAY_DIR"] #get the delaydir, what I actually want to use for updating phase_dir
    newphasedir = delaydir
    #This means that can run this code on a MS multiple times and it will work
    string_delay = format_phasedir(delaydir)
    string_phase = format_phasedir(newphasedir)
    print 'The delay direction is {0} and the new phase direction is {1}'.format(string_delay,string_phase)
    print ''
    #print to check against delay direction and make sure is the same
    
    #do the actual update
    t_field.putcell("PHASE_DIR", 0, newphasedir)  #update the phase direction
    t_field.flush() #make sure changes are written to MS
    
    #recalculate uv coordinates
    pt.taql('update {0} set UVW = mscal.newuvw()'.format(msfile)) 
    