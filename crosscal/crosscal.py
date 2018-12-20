#python "module" for tests and inventories related to
#confirming the ability and calculating appropriate scaling
#for transfering solutions from central beam to other beams

#load necessary packages
import os
import numpy as np
from astropy.io import ascii
import sys
sys.path.append('/home/adams/commissioning/general')
import apertif_utility as aputil
import apercal
import casacore.tables as pt
import matplotlib.pyplot as plt

#define a class that is used for specifying the scans to be investigated

"""CLASS definitions"""

class ScanSpecification(object):
    def __init__(self):
        #set all attributes to empty strings to initialize
        #that way I can check for what is set later
        self.startdate=''
        self.enddate=''
        self.beam=''
        self.startscan=''
        self.endscan=''
        self.nscan=''
        self.badscans = ''
    def setstartdate(self,sd):
        #want to do some checking that startdate is in the right format
        #turn into a string
        sdstr=str(sd)
        if len(sdstr) != 8:
            print "Startdate must be  of format 'YYYYMMDD'"
        elif sdstr[0:3] != '201':
            print "Startdate must be a string of format 'YYYYMMDD', starting from 2010"
        else:
            self.startdate=sdstr
    def setenddate(self,ed):
        edstr=str(ed)
        if len(edstr) != 8:
            print "Enddate must be of format 'YYYYMMDD'"
        elif edstr[0:3] != '201':
            print "Enddate must be of format 'YYYYMMDD', starting from 2010"
        else:
            self.enddate = edstr
    def setbeam(self,bm):
        bmstr=str(bm)
        if len(bmstr) == 1:
            newbmstr = '{0}{1}'.format(0,bmstr)
            self.beam = newbmstr
        elif float(bmstr) != int(float(bmstr)):
            print "Beam must be an integer"
        elif len(bmstr) != 2:
            print "Beam must be of format 'NN'"
        elif float(bmstr) >39 or float(bmstr) <0:
            print "Beam must be an integer between 0 and 39"
        else:
            self.beam = bmstr
    def setstartscan(self,stsc):
        scstr = str(stsc)
        #Assume format is YYMMDDXXX 
        #this is valid from Feb 2018
        #Jan and earlier had only two scan numbers
        #I think just focusing on that date range (at least to start) is okay
        if len(scstr) != 9:
            print "Start scan must be of format 'YYMMDDXXX'"
        elif scstr[0:2] != '18': #can add more years later if need be
            print "Start scan must be of format 'YYMMDDXXX', with a year of 18"
        else:
            self.startscan = scstr
    def setendscan(self,edsc):
        scstr = str(edsc)
        #Assume format is YYMMDDXXX 
        #this is valid from Feb 2018
        #Jan and earlier had only two scan numbers
        #I think just focusing on that date range (at least to start) is okay
        if len(scstr) != 9:
            print "Start scan must be of format 'YYMMDDXXX'"
        elif scstr[0:2] != '18': #can add more years later if need be
            print "Start scan must be of format 'YYMMDDXXX', with a year of 18"
        else:
            self.endscan = scstr    
    def setnscan(self,ns):
        nsstr = str(ns)
        #want to check if it is a number
        #but doesn't work how I expect
        try:
            float(nsstr)
            if float(nsstr) == int(float(nsstr)):
                self.nscan = nsstr
            else:
                new_ns = str(int(float(nsstr)))
                self.nscan = new_ns
        except ValueError:
            print 'nscan must be a number'
    def setbadscans(self,badscanlist):
        if type(badscanlist) is list:
            newbadscanlist = []
            for badsc in badscanlist:
                #check each element of list
                badscstr = str(badsc)
                if len(badscstr) != 9:
                    print "Bad scan must be of format 'YYMMDDXXX'"
                elif badscstr[0:2] != '18': #can add more years later if need be
                    print "Bad scan must be of format 'YYMMDDXXX', with a year of 18"
                else:
                    newbadscanlist.append(badscstr)
            self.badscans = newbadscanlist  #new list is strings and only has properly formatted scans                
        else:
            print "Bad scans must be provided as a list"

            

"""Prepare data"""

def get_scan_list(scans,obsrecordfile):
    #This takes a ScanObject that specifies the scans wanted
    #Use this information to produce a scan and beam list
    #All data is on ALTA now, so need to search there
    #ATDB is not yet up to the task, so use Apertif Observation Record
    #provide this as a specific input in csv format
    
    #provide a list of bad scans for failed data
    #these will be skipped in creating lists
    
    #first, read the observation record in:
    print 'Reminder: Make sure observation record is up-to-date!' #reminder
    obsrecord = ascii.read(obsrecordfile)

    #start by checking for different modes of operation
    #two modes will be 
    #(1) switch: for beam switchign observations in a row
    #(2) variability: for tracking long-term variability in a given scan
    if scans.startdate != '':
        if scans.enddate != '':
            if scans.beam != '':
                mode = 'variability'
                print 'In variability scan mode'
            else:
                print 'Variability mode is missing a beam specification'
                mode=None
        else:
            print 'Variability mode is missing an enddate'
            mode = None
    elif scans.startscan !='':
        if scans.endscan !='':
            mode='switch'
            print 'In switching scan mode'
        elif scans.nscan !='':
            mode='switch'
            print 'In switching mode'
        else:
            print 'Switching mode needs an end scan or nscan specification'
            mode = None
    else:
        mode = None
        print "No scan specification set"
    
    #set placeholders for scan and beam list        
    scan_list = []
    beam_list = []


    if mode == 'switch':
        if scans.nscan == '':
            #if nscan isn't set, calculate it
            nscans = int(scans.endscan) - int(scans.startscan) + 1
            scans.setnscan(nscans)
            #find number of bad scans and account for that
        scan_list = np.empty(int(scans.nscan),dtype=object) #set scanlist to be lenght of nscan
        beam_list = np.empty(int(scans.nscan),dtype=object) #set beamlist to length of nscan
        #iterate through each scan and add it to the scan_list
        #parse obsrecord to find relevant beam aand add beam to beamlist
        badind = [] #track indices for bad scans
        for n in xrange(int(scans.nscan)):
            scannumber = int(scans.startscan)+n #get scan number as an int
            #check if bad, if yes, mark it:
            if str(scannumber) in scans.badscans:
                badind.append(n)
            scan_list[n] = str(scannumber)  #write everything as strings
            #find the entry in obsrecord TaskId and parse beam from Data entry
            ind = np.where(obsrecord['TaskId'] == scannumber)[0]
            name = str(obsrecord['Data'][ind]) 
            #have to turn into a string because this is 'MaskedColumn'
            name_split = name.split('_')
            #find the last entry for beam number.
            #There should always be three entries - scan (with column name, annoyingly), source and beam
            #If there are only two entries, assume beam is 0. 
            #This might not have been true in past but should be always now.
            if len(name_split) <3:
                beam = 0
            else:
                beam = int(name_split[-1])
            beam_list[n] = str(beam)  #write as string
            
    if mode == 'variability':
        #need to implement this but it isn't critical yet.
        print "Identifying scanlist and beamlist for variability mode isn't implemented yet"

    #delete bad scans from list
    #not the most elegant way but works for now
    if len(badind) >0:
        scans = np.delete(scan_list,badind)
        beams = np.delete(beam_list,badind)
    else:
        scans = scan_list
        beams = beam_list

    return mode,scans,beams



def copy_scans(scans,obsrecordfile,basedir,run=False):
    #This will use the output of get_scan_list 
    #to retrieve data from ALTA and move to targetdest
    #I will take same input as get_scan_list and run it here
    #Want to minimize number of calls by user
    #Maybe the proper way to do this is to have global variables, 
    #But that's a step too complicated for me, I think
    mode,scan_list,beam_list = get_scan_list(scans,obsrecordfile)
    #get the values I need
    
    #get the data
    

    if run==False:
        print 'In verification mode, will print commands to screen'
    elif run==True:
        print 'In running mode, will execute commands'
    else:
        print 'Mode is not defined. Choose verify or run'
    
    #now iterate through scan/beam lists and run command to retrieve data
    
    for scan,beam in zip(scan_list,beam_list):
        #first, move into targetdest (create if needed)
        #need a separate directory for each dataset
        #bleh
        #name by scan number (best for when I add other mode)
        targetdir = '{0}/{1}/00/raw'.format(basedir,scan)
        if os.path.exists(targetdir):
            pass
        else:
            os.makedirs(targetdir)
        #then change to working directory to copy data
        #think this is how the usage of ALTA data module from Vanessa is
        os.chdir(targetdir)
        print 'Moved to and copying data to {0}'.format(os.getcwd())
    
        #first parse scan as required:
        scandate = scan[0:6]
        scannum = scan[6:9]
        #format string command to match usage:
            # ALTA data transfer: Uses the iROD client to transfer data from ALTA
            # Example usage: >> python getdata_alta.py 180316 004-010 00-36
            # V.A. Moss (vmoss.astro@gmail.com)
        #hope is that since I did sys.path.append('/home/adams/altadata'),
        #I don't have to specify full path
        string_command = "python /home/adams/altadata/getdata_alta.py {0} {1}-{1} {2:0>2}-{2:0>2}".format(scandate,scannum,beam)
        #!!!!!!!Important question - do I have specify beams with double integers?
        #Thayt's annoying because it doesn't match naming scheme
        #But maybe I can do with format - yep!
        #check if data already exists - if it does, don't do anything (set run=False to print to screen)
        filename = "WSRTA{0}_B{1:0>3}.MS".format(scan,beam)
        filerun = True
        if os.path.exists(filename):
            filerun=False #file exists, so don't need to copy again            
        if run==False:
            print string_command
        if filerun == False:
            print 'File already exists'
            #reset filerun
            filerun=True
        elif run==True:
            #and run the command:
            os.system(string_command)
            print 'Running '+string_command

def fix_source_name(scans,obsrecordfile,basedir,flip=False,unflip=False):
    #Source names have beam in them (critical for get_scans!)
    #But that will cause calibration to crash
    #CASA isn't very smart
    #So I'll have to be smart instead and update everything
    #I'll also give option here of flipping coordinates
    #Since I'm already messing w/ MS, this seems appropriate place to do so.
    #and will also include unflip flag, in case I mess things up
    mode,scan_list,beam_list = get_scan_list(scans,obsrecordfile)
    for scan,beam in zip(scan_list,beam_list):
        targetdir = '{0}/{1}/00/raw'.format(basedir,scan)
        msfile = "{0}/WSRTA{1}_B{2:0>3}.MS".format(targetdir,scan,beam)
        t_field = pt.table(msfile+"::FIELD", readonly=False)
        name = t_field[0]['NAME']
        name_split = name.split('_') #split by underscore - works also if no underscore
        t_field.putcell("NAME", 0, name_split[0])  #update source name to anything before first underscore (or original name if no underscore)
        t_field.flush()
        if (flip == True) and (unflip == True):
            print "Both flipping and unflipping. No net change"
        if flip == True:
            aputil.flip_ra(msfile)
        if unflip == True:
            aputil.unflip_ra(msfile)
    
    
def flag_scans(scans,obsrecordfile,basedir,cfgfile,edges=True,ghosts=True):
    #This will use output of get_scan_list
    #plus apercal.preflag() to flag the data
    #basedir should be in cfg file but specify manually to make sure things match
    mode,scan_list,beam_list = get_scan_list(scans,obsrecordfile)
    preflag = apercal.preflag(cfgfile)
    preflag.preflag_edges = edges #option to turn on/off
    preflag.preflag_ghosts = ghosts #option to turn on/off
    #set source and pol cal to False - will iterate through updating fluxcal and running preflag one at a time
    preflag.preflag_manualflag_polcal = False
    preflag.preflag_manualflag_target = False
    preflag.preflag_aoflagger_polcal = False 
    preflag.preflag_aoflagger_target = False
    #should only find central beam, but specify anyway
    preflag.preflag_aoflagger_targetbeams = '00' 
    #now iterate through all my scans/beams:
    for scan,beam in zip(scan_list,beam_list):
        preflag.fluxcal = "WSRTA{0}_B{1:0>3}.MS".format(scan,beam)
        preflag.basedir = "{0}/{1}/".format(basedir,scan) #basedir for Apercal is different than my basedir, needs trailing slash
        print "Setting fluxcal to WSRTA{0}_B{1:0>3}.MS".format(scan,beam)
        #treat everything as fluxcal - this is for purpose of iterating and flagging
        #may have to be more careful when deriving solutions - will care about flux vs. pol eventually then
        #and run preflag
        print "Flagging data set {2}00/raw/WSRTA{0}_B{1:0>3}.MS".format(scan,beam,preflag.basedir)
        preflag.go()


    
       
def calibrate_scans(scans,obsrecordfile,basedir,cfgfile):
    #This will use scan_list,beam_list from get_scan_list
    #plus apercal.ccal to calibrate the data
    #HAVE TO CHANGE SOURCE NAME! - this is done in fix_source_name
    print "Warning! Have you run fix_source_name?"
    mode,scan_list,beam_list = get_scan_list(scans,obsrecordfile)
    
    #load ccal module
    ccal = apercal.ccal(cfgfile)
    ccal.crosscal_transfer_to_target = False #there is no target
    for scan,beam in zip(scan_list,beam_list):
        ccal.fluxcal = "WSRTA{0}_B{1:0>3}.MS".format(scan,beam)
        ccal.basedir = "{0}/{1}/".format(basedir,scan) #basedir for Apercal is different than my basedir, plus trailing slash
        print "Setting fluxcal to WSRTA{0}_B{1:0>3}.MS".format(scan,beam)
        print "Calibrating data set {2}00/raw/WSRTA{0}_B{1:0>3}.MS".format(scan,beam,ccal.basedir)
        ccal.go()

def clear_calibration_scans(scans,obsrecordfile,basedir,cfgfile):  
    #this will clear calibration solution
    #that way can rerun w/out having to redo everything
    mode,scan_list,beam_list = get_scan_list(scans,obsrecordfile)
    
    #load ccal module
    ccal = apercal.ccal(cfgfile)
    for scan,beam in zip(scan_list,beam_list):
        ccal.fluxcal = "WSRTA{0}_B{1:0>3}.MS".format(scan,beam)
        ccal.basedir = "{0}/{1}/".format(basedir,scan)
        ccal.clear()

"""Bandpass solutions"""    
    
def get_bp_sols(bptable):
    #takes a path to a BP solution table
    #This should also work for the delay table - same keyword and same structure (one row per ant)
    #Gain tables don't work - values per time
    #K tables don't work - FPARAM keyword
    taql_command = ("SELECT TIME,abs(CPARAM) AS amp, arg(CPARAM) AS phase, "
                    "FLAG FROM {0}").format(bptable)
    t=pt.taql(taql_command)
    times = t.getcol('TIME')
    amp_sols=t.getcol('amp')
    phase_sols = t.getcol('phase')
    flags = t.getcol('FLAG')

    taql_antnames = "SELECT NAME FROM {0}::ANTENNA".format(bptable)
    t= pt.taql(taql_antnames)
    ant_names=t.getcol("NAME")

    taql_freq = "SELECT CHAN_FREQ FROM {0}::SPECTRAL_WINDOW".format(bptable)
    t = pt.taql(taql_freq)
    freqs = t.getcol('CHAN_FREQ')

    return ant_names,times,freqs,amp_sols,phase_sols,flags

def compare_scan_solution_bp(scans,obsrecordfile,basedir,norm=True,refscan=''):
    #This will compare scan BP solutions to reference (given) for amp & phase 
    #(easiest to do together at once)
    #will return a multi-d array (ant, freq, pol, scan)
    #If norm=True, amp solution is divided by reference and phase ref is subtracted from sol
    #otherwise it returns all the solutions into array
    
    #first, get scan and beam list:
    mode,scan_list,beam_list = get_scan_list(scans,obsrecordfile)
    
    #check if want normalized solutions
    #if so, get reference values
    if norm == True:        
        if refscan =='':
            print 'Must provide a reference scan for normalization!'
            print 'Setting normalization to False'
            norm = False #set normalization to False since there is no ref scan
            ref_amp_sol = np.nan
            ref_phase_sol = np.nan
        else:
            print 'Will normalize solutions'
            ind = np.where(str(refscan) == scan_list)[0]
            scan = scan_list[ind][0]
            beam = beam_list[ind][0]
            refbpsol = "{0}/{1}/00/raw/WSRTA{1}_B{2:0>3}.Bscan".format(basedir,scan,beam)
            #assumes naming convention in Apercal won't change!
            ant_names,times,freqs,ref_amp_sol,ref_phase_sol,flags = get_bp_sols(refbpsol)
    
    #now iterate through each beam
    #will want to normalize by reference, if that option is set
    #create array to hold everything before I start
    #easy if I have reference, more difficult otherwise
    #So maybe just get first value as a test
    testsol = "{0}/{1}/00/raw/WSRTA{1}_B{2:0>3}.Bscan".format(basedir,scan_list[0],beam_list[0])
    ant_names,times,freqs,amps,phases,flags = get_bp_sols(testsol)
    bp_amp_vals = np.empty((amps.shape[0],amps.shape[1],amps.shape[2],int(scans.nscan)))
    bp_phase_vals = np.empty((phases.shape[0],phases.shape[1],phases.shape[2],int(scans.nscan)))
    #iterate through scans:
    for n,(scan,beam) in enumerate(zip(scan_list,beam_list)):
        bpsol = "{0}/{1}/00/raw/WSRTA{1}_B{2:0>3}.Bscan".format(basedir,scan,beam)
        ant_names,times,freqs,amps,phases,flags = get_bp_sols(bpsol)
        if norm == True:
            bp_amp = amps / ref_amp_sol
            bp_phase = phases - ref_phase_sol 
            #"nromalize" phase by subtracting - care about absolute deviation
        else:
            bp_amp = amps
            bp_phase = phases
        #check for flags and mask
        bp_amp[flags] = np.nan
        bp_phase[flags] = np.nan
        #populate arrays
        bp_amp_vals[:,:,:,n] = bp_amp
        bp_phase_vals[:,:,:,n] = bp_phase * 180/np.pi #put in degrees
               
    return ant_names,times,freqs,bp_amp_vals,bp_phase_vals


def plot_compare_bp_beam(scans,obsrecordfile,basedir,norm=True,
                         refscan='',plotmode='amp',pol=0,nx=3,ymin=0,ymax=0,plotsize=4,
                        figname=''):
    #this will generate plots that compare BP solutions between beams
    #It can run with option amplitude or phase, default amp
    #Defaults to showing pol 0, can also change
    
    mode,scan_list,beam_list = get_scan_list(scans,obsrecordfile)
    (ant_names,times,freqs,
     bp_amp_vals,bp_phase_vals) = compare_scan_solution_bp(scans,obsrecordfile,basedir,
                                                           norm=norm,refscan=refscan)        
    if plotmode == 'amp':
        plotvals = bp_amp_vals
        print 'plotting amplitude solutions'
    elif plotmode == 'phase':
        plotvals = bp_phase_vals
        print 'plotting phase solutions'
    else:
        print 'plotmode={} not defined'.format(mode)
    
    #now setup the plotting environment
    #total number of plots is number of scans
    nplots = len(scan_list)
    ny = int(np.ceil(nplots/float(nx)))
    #and set up size that I will want
    #say 3 inches per plot
    xsize = nx*plotsize
    ysize = ny*plotsize
    #and I want global limits, for best comparison
    #this could potentially change
    #will have a total offset of 
    if ymin == 0:
        ymin = np.nanmin(plotvals)
    if ymax == 0:
        ymax = np.nanmax(plotvals)
    xmin = np.nanmin(freqs)
    xmax = np.nanmax(freqs)
    
    #create figure object
    #this is what I return from program
    fig= plt.figure(figsize=(xsize,ysize))
    plt.xlim(xmin,xmax) # Limit the plot to the minimum and maximum frequencies
    plt.suptitle('Bandpass {0}, normalization {1}'.format(plotmode,norm), fontsize='large')


    for n,(scan,beam) in enumerate(zip(scan_list,beam_list)):
        plt.subplot(ny, nx, n+1)
        for a,ant in enumerate(ant_names):
            plt.scatter(freqs[0,:], plotvals[a,:,0,n],label=ant,
                       marker=',',s=0.5)
#        for sol in range(plotvals.shape[2]):
#            plt.plot(plotfreqs, (plotvals[a,:,sol] + np.full(plotvals.shape[1],offset*sol)))
        plt.title('{0}, beam {1}'.format(scan,beam))
        plt.xlim(xmin,xmax) # Limit the plot to the minimum and maximum frequencies
        plt.ylim(ymin,ymax)
           
    plt.legend()
    
    if figname != '':
        plt.savefig(figname)
        
    return fig


"""Gain solutions"""

def get_gain_sols(gaintable):
    #takes path to an gain solution table
    #returns amp,phase as function of time (no frequency info)
    
    #first get antenna names, need to iterate over later
    taql_antnames = "SELECT NAME FROM {0}::ANTENNA".format(gaintable)
    t= pt.taql(taql_antnames)
    ant_names=t.getcol("NAME")
    
    #then get number of times
    #need this for setting shape
    taql_time =  "select TIME from {0} orderby unique TIME".format(gaintable)
    t= pt.taql(taql_time)
    times = t.getcol('TIME') 
    
    #then iterate over antenna
    #set array sahpe to be [n_ant,n_time,n_stokes]
    #how can I get n_stokes? Could be 2 or 4, want to find from data
    #get 1 data entry
    taql_stokes = "SELECT abs(CPARAM) AS amp from {0} limit 1" .format(gaintable)
    t_pol = pt.taql(taql_stokes)
    pol_array = t_pol.getcol('amp')
    n_stokes = pol_array.shape[2] #shape is time, one, nstokes
    
    amp_ant_array = np.empty((len(ant_names),len(times),n_stokes),dtype=object)
    phase_ant_array = np.empty((len(ant_names),len(times),n_stokes),dtype=object)
    flags_ant_array = np.empty((len(ant_names),len(times),n_stokes),dtype=bool)
    
    for ant in xrange(len(ant_names)):
        taql_command = ("SELECT abs(CPARAM) AS amp, arg(CPARAM) AS phase, FLAG FROM {0} " 
                        "WHERE ANTENNA1={1}").format(gaintable,ant)
        t = pt.taql(taql_command)
        amp_ant_array[ant,:,:] = t.getcol('amp')[:,0,:]
        phase_ant_array[ant,:,:] = t.getcol('phase')[:,0,:]
        flags_ant_array[ant,:,:] = t.getcol('FLAG')[:,0,:]

    return ant_names,times,amp_ant_array,phase_ant_array,flags_ant_array


def compare_scan_solution_gain(scans,obsrecordfile,basedir,norm=True,refscan=''):
    #This will collect all gain solutions for a scan of beams
    #If set, will nromalize to a reference.
    #Note that each scan will have different times, 
    #so will compute a scalar average
    #and use that for comparison
    #for now, assume separated by beam but by time should also work
    
    #first, get scan and beam list:
    mode,scan_list,beam_list = get_scan_list(scans,obsrecordfile)
    
    #check if want normalized solutions
    #if so, get reference values
    if norm == True:        
        if refscan =='':
            print 'Must provide a reference scan for normalization!'
            print 'Setting normalization to False'
            norm = False #set normalization to False since there is no ref scan
            ref_amp_sol = np.nan
            ref_phase_sol = np.nan
        else:
            print 'Will normalize solutions'
            ind = np.where(str(refscan) == scan_list)[0]
            scan = scan_list[ind][0]
            beam = beam_list[ind][0]
            refsol = "{0}/{1}/00/raw/WSRTA{1}_B{2:0>3}.G1ap".format(
                basedir,scan,beam)
            #assumes naming convention in Apercal won't change!
            (ant_names,times,amp_ant_array,phase_ant_array,
             flag_ant_array) = get_gain_sols(refsol)
            #average over time
            #have array of shape (nant,nstokes):
            refamp_ant_array = np.mean(amp_ant_array,axis=1)
            refphase_ant_array = np.mean(phase_ant_array,axis=1)
            array_shape = refamp_ant_array.shape
            refamp = refamp_ant_array.reshape(array_shape[0],1,array_shape[1])
            refphase = refphase_ant_array.reshape(array_shape[0],1,array_shape[1])
    
    #now iterate through each beam
    #will want to normalize by reference, if that option is set
    #create array to hold everything before I start
    #easy if I have reference, more difficult otherwise
    #So maybe just get first value as a test
    """
    #Note that I will set length of times based on this, 
    #But this may not be true for all scans
    #So I will force everything to same shape,
    #Or I could try using dictionaries
    #or objects instead
    #that's probably the best way to do it
    #but i'll go easy/brute-force for now
    #and leave that for future improvement
    """
    testsol = "{0}/{1}/00/raw/WSRTA{1}_B{2:0>3}.G1ap".format(
        basedir,scan_list[0],beam_list[0])
    ant_names,times,amp_sols,phase_sols,flags = get_gain_sols(testsol)
    gain_amp_vals = np.empty((amp_sols.shape[0],amp_sols.shape[1],
                              amp_sols.shape[2],int(scans.nscan)))
    gain_phase_vals = np.empty((phase_sols.shape[0],phase_sols.shape[1],
                                phase_sols.shape[2],int(scans.nscan)))
    time_vals = np.empty((amp_sols.shape[1],int(scans.nscan)))
    ntimes = amp_sols.shape[1]
    #iterate through scans:
    for n,(scan,beam) in enumerate(zip(scan_list,beam_list)):
        gainsol = "{0}/{1}/00/raw/WSRTA{1}_B{2:0>3}.G1ap".format(
            basedir,scan,beam)
        ant_names,times,amp_sols,phase_sols,flags = get_gain_sols(gainsol)
        if norm == True:
            #check for length of time axis
            if len(times) == ntimes:
                gain_amp = np.divide(amp_sols, refamp)
                gain_phase = np.subtract(phase_sols,refphase)
                gain_times = times
                gain_flags = flags
                #"nromalize" phase by subtracting - care about absolute deviation
            elif len(times) > ntimes:
                gain_amp = np.divide(amp_sols[:,0:ntimes,:], refamp)
                gain_phase = np.subtract(phase_sols[:,0:ntimes,:],refphase)
                gain_times = times[0:ntimes]
                gain_flags = flags[0:ntimes]
                #"nromalize" phase by subtracting - care about absolute deviation
            else:
                gain_amp = np.full((amp_sols.shape[0],ntimes,amp_sols.shape[2]),
                                   np.nan)
                gain_phase = np.full((amp_sols.shape[0],ntimes,amp_sols.shape[2]),
                                     np.nan)
                gain_flags = np.full((amp_sols.shape[0],ntimes,amp_sols.shape[2]),
                                     True)
                gain_times = np.full(ntimes,np.nan)
                gain_amp[:,0:len(times),:] = np.divide(amp_sols, refamp)
                gain_phase[:,0:len(times),:] = np.subtract(phase_sols,refphase)
                gain_times[0:len(times)] = times
                gain_flags[:,0:len(times),:] = flags
        else:
            #check for length of time axis:
            if len(times) == ntimes:
                gain_amp = amp_sols
                gain_phase = phase_sols
                gain_times = times
                gain_flags = flags
            elif len(times) > ntimes:
                gain_amp = amp_sols[:,0:ntimes,:]
                gain_phase = phase_sols[:,0:ntimes,:]
                gain_times = times[0:ntimes]
                gain_flags = flags[0:ntimes]
            else:
                gain_amp = np.full((amp_sols.shape[0],ntimes,amp_sols.shape[2]),
                                   np.nan)
                gain_phase = np.full((amp_sols.shape[0],ntimes,amp_sols.shape[2]),
                                     np.nan)
                gain_flags = np.full((amp_sols.shape[0],ntimes,amp_sols.shape[2]),
                                      True)
                gain_times = np.full(ntimes,np.nan)
                gain_amp[:,0:len(times),:] = amp_sols
                gain_phase[:,0:len(times),:] = phase_sols
                gain_times[0:len(times)] = times
                gain_flags[:,0:len(times),:] = flags

        #find and replace flags
        #print scan, gain_amp.shape,gain_flags.shape
        gain_amp[gain_flags] = np.nan
        gain_phase[gain_flags] = np.nan
        #gain_times[gain_flags] = np.nan
        #populate arrays
        gain_amp_vals[:,:,:,n] = gain_amp
        gain_phase_vals[:,:,:,n] = gain_phase * 180/np.pi #put in degrees
        time_vals[:,n] = gain_times
        
    return ant_names,time_vals,gain_amp_vals,gain_phase_vals


def plot_compare_gain_beam(scans,obsrecordfile,basedir,
                           norm=True,refscan='',plotmode='amp',
                           pol=0,nx=3,ymin=0,ymax=0,plotsize=4,
                          figname=''):
    #this will generate plots that compare BP solutions between beams
    #It can run with option amplitude or phase, default amp
    #Defaults to showing pol 0, can also change
    
    mode,scan_list,beam_list = get_scan_list(scans,obsrecordfile)
    (ant_names,time_vals,gain_amp_vals,
     gain_phase_vals) = compare_scan_solution_gain(scans,
                                                   obsrecordfile,
                                                   basedir,norm=norm,
                                                   refscan=refscan)

    if plotmode == 'amp':
        plotvals = gain_amp_vals
        print 'plotting amplitude solutions'
    elif plotmode == 'phase':
        plotvals = gain_phase_vals
        print 'plotting phase solutions'
    else:
        print 'plotmode={} not defined'.format(mode)
    
    #now setup the plotting environment
    #total number of plots is number of scans
    nplots = len(scan_list)
    ny = int(np.ceil(nplots/float(nx)))
    #and set up size that I will want
    #say 3 inches per plot
    xsize = nx*plotsize
    ysize = ny*plotsize
    #and I want global limits, for best comparison
    #this could potentially change
    #will have a total offset of 
    if ymin == 0:
        ymin = np.nanmin(plotvals)
    if ymax == 0:
        ymax = np.nanmax(plotvals)
    
    #create figure object
    #this is what I return from program
    fig= plt.figure(figsize=(xsize,ysize))
    plt.suptitle('Gain {0}, normalization {1}'.
                 format(plotmode,norm), fontsize='large')


    for n,(scan,beam) in enumerate(zip(scan_list,beam_list)):
        plt.subplot(ny, nx, n+1)
        for a,ant in enumerate(ant_names):
            plt.plot(time_vals[:,n], plotvals[a,:,pol,n],label=ant)
#        for sol in range(plotvals.shape[2]):
#            plt.plot(plotfreqs, (plotvals[a,:,sol] + np.full(plotvals.shape[1],offset*sol)))
        plt.title('{0}, beam {1}'.format(scan,beam))
        plt.ylim(ymin,ymax)
           
    plt.legend()
    
    if figname != '':
        plt.savefig(figname)
    
    return fig


"""Model data"""

def get_model(msfile):
    #take a MS file and get the model column
    #average over baselines/times
    #want amp/phase vs. frequency
    taql_command = "SELECT abs(gmeans(MODEL_DATA)) AS amp, arg(gmeans(MODEL_DATA)) AS phase FROM {0}".format(msfile)
    t = pt.taql(taql_command)
    amp = t.getcol('amp')[0,:,:]
    phase = t.getcol('phase')[0,:,:]
    taql_freq = "SELECT CHAN_FREQ FROM {0}::SPECTRAL_WINDOW".format(msfile)
    t = pt.taql(taql_freq)
    freqs = t.getcol('CHAN_FREQ')[0,:]
    
    return freqs,amp,phase

def compare_scan_model(scans,obsrecordfile,basedir):
    #take a scan object and get the model for each observation
    #this will confirm that sources are properly named
    
    #first, get scan and beam list:
    mode,scan_list,beam_list = get_scan_list(scans,obsrecordfile)
    
    #now iterate through each beam
    #start with a test to get dimensions
    
    testscan = "{0}/{1}/00/raw/WSRTA{1}_B{2:0>3}.MS".format(
        basedir,scan_list[0],beam_list[0])
    freqs,amp,phase = get_model(testscan)
    model_amp_array = np.empty((amp.shape[0],amp.shape[1],int(scans.nscan)))
    model_phase_array = np.empty((phase.shape[0],phase.shape[1],int(scans.nscan)))

    #iterate through scans:
    for n,(scan,beam) in enumerate(zip(scan_list,beam_list)):
        scan = "{0}/{1}/00/raw/WSRTA{1}_B{2:0>3}.MS".format(
            basedir,scan,beam)
        freqs,amp,phase = get_model(scan)
        model_amp_array[:,:,n] = amp
        model_phase_array[:,:,n] = phase
        
    return freqs,model_amp_array,model_phase_array

def plot_compare_scan_model(scans,obsrecordfile,basedir,plotmode='amp',
                            pol=0,nx=3,ymin=0,ymax=0,plotsize=4,
                           figname=''):
    #plot model for each scan
    mode,scan_list,beam_list = get_scan_list(scans,obsrecordfile)
    (freqs,model_amp_array,
     model_phase_array) = compare_scan_model(scans,obsrecordfile,
                                             basedir)
    #figure out plotmode
    if plotmode == 'amp':
        plotvals = model_amp_array
        print 'plotting model amplitudes'
    elif plotmode == 'phase':
        plotvals = model_phase_array
        print 'plotting model phases'
    else:
        print 'plotmode={} not defined'.format(plotmode)
        
    #now setup the plotting environment
    #total number of plots is number of scans
    nplots = len(scan_list)
    ny = int(np.ceil(nplots/float(nx)))
    xsize = nx*plotsize
    ysize = ny*plotsize

    if ymin == 0:
        ymin = np.nanmin(plotvals)
    if ymax == 0:
        ymax = np.nanmax(plotvals)
    xmin = np.nanmin(freqs)
    xmax = np.nanmax(freqs)
    
    #create figure object
    #this is what I return from program
    fig= plt.figure(figsize=(xsize,ysize))
    plt.xlim(xmin,xmax) # Limit the plot to the minimum and maximum frequencies
    plt.suptitle('Model {0}'.format(plotmode), fontsize='large')

    for n,(scan,beam) in enumerate(zip(scan_list,beam_list)):
        plt.subplot(ny, nx, n+1)
        plt.scatter(freqs,plotvals[:,pol,n],marker=',',s=0.5)
        plt.title('{0}, beam {1}'.format(scan,beam))
        plt.xlim(xmin,xmax) # Limit the plot to the minimum and maximum frequencies
        plt.ylim(ymin,ymax)
        
    if figname != '':
        plt.savefig(figname)
        
    return fig
    
    
"""Calibrated data"""

def get_calibrated_data(msfile):
    #first get antenna names, need to iterate over later
    taql_antnames = "SELECT NAME FROM {0}::ANTENNA".format(msfile)
    t= pt.taql(taql_antnames)
    ant_names=t.getcol("NAME")

    #then get frequencies:
    taql_freq = "SELECT CHAN_FREQ FROM {0}::SPECTRAL_WINDOW".format(msfile)
    t = pt.taql(taql_freq)
    freqs = t.getcol('CHAN_FREQ')[0,:]
    
    #and number of stokes params
    taql_stokes = "SELECT abs(DATA) AS amp from {0} limit 1" .format(msfile)
    t_pol = pt.taql(taql_stokes)
    pol_array = t_pol.getcol('amp')
    n_stokes = pol_array.shape[2] #shape is time, one, nstokes
    
    #take MS file and get calibrated data
    amp_ant_array = np.empty((len(ant_names),len(freqs),n_stokes),dtype=object)
    phase_ant_array = np.empty((len(ant_names),len(freqs),n_stokes),dtype=object)
    
    for ant in xrange(len(ant_names)):
        taql_command = ("SELECT abs(gmeans(CORRECTED_DATA[FLAG])) AS amp, "
                        #"gmeans(FLAG), "
                        "arg(gmeans(CORRECTED_DATA[FLAG])) AS phase FROM {0} "
                        "WHERE ANTENNA1!=ANTENNA2 && "
                        "(ANTENNA1={1} || ANTENNA2={1})").format(msfile,ant)
        #use FLAG as mask to only use unflagged data
        t = pt.taql(taql_command)
        test=t.getcol('amp')
        amp_ant_array[ant,:,:] = t.getcol('amp')[0,:,:]
        phase_ant_array[ant,:,:] = t.getcol('phase')[0,:,:]
    
    return ant_names,freqs,amp_ant_array,phase_ant_array


def compare_scan_calibrated_data(scans,obsrecordfile,basedir,norm=True,refscan=''):
    #can normalize to reference
    #useful to see if there are systematic issues, but wouldn't expect taht
    
    #first, get scan and beam list:
    mode,scan_list,beam_list = get_scan_list(scans,obsrecordfile)
    
    #check if want normalized solutions
    #if so, get reference values
    if norm == True:        
        if refscan =='':
            print 'Must provide a reference scan for normalization!'
            print 'Setting normalization to False'
            norm = False #set normalization to False since there is no ref scan
            ref_amp = np.nan
            ref_phase = np.nan
        else:
            print 'Will normalize solutions'
            ind = np.where(str(refscan) == scan_list)[0]
            scan = scan_list[ind][0]
            beam = beam_list[ind][0]
            ref = "{0}/{1}/00/raw/WSRTA{1}_B{2:0>3}.MS".format(basedir,scan,beam)
            #assumes naming convention in Apercal won't change!
            ant_names,freqs,ref_amp,ref_phase = get_calibrated_data(ref)

    #now iterate through each beam
    #first set dimenstions
    testsol = "{0}/{1}/00/raw/WSRTA{1}_B{2:0>3}.MS".format(basedir,scan_list[0],beam_list[0])
    ant_names,freqs,amps,phases = get_calibrated_data(testsol)
    amp_vals = np.empty((amps.shape[0],amps.shape[1],amps.shape[2],int(scans.nscan)))
    phase_vals = np.empty((phases.shape[0],phases.shape[1],phases.shape[2],int(scans.nscan)))
    #iterate through scans:
    for n,(scan,beam) in enumerate(zip(scan_list,beam_list)):
        data = "{0}/{1}/00/raw/WSRTA{1}_B{2:0>3}.MS".format(basedir,scan,beam)
        ant_names,freqs,amps,phases = get_calibrated_data(data)
        if norm == True:
            data_amp = amps / ref_amp
            data_phase = phases - ref_phase #"nromalize" phase by subtracting - care about absolute deviation
        else:
            data_amp = amps
            data_phase = phases

        amp_vals[:,:,:,n] = data_amp
        phase_vals[:,:,:,n] = data_phase * 180/np.pi #put in degrees
        
    return ant_names,freqs,amp_vals,phase_vals        
            

def plot_compare_calibrated_data_beam(scans,obsrecordfile,basedir,norm=True,
                                      refscan='',plotmode='amp',pol=0,nx=3,ymin=0,ymax=0,plotsize=4,
                                     figname=''):
    #this will generate plots that compare BP solutions between beams
    #It can run with option amplitude or phase, default amp
    #Defaults to showing pol 0, can also change
    
    mode,scan_list,beam_list = get_scan_list(scans,obsrecordfile)
    ant_names,freqs,amp_vals,phase_vals = compare_scan_calibrated_data(scans,obsrecordfile,basedir,
                                                                       norm=norm,refscan=refscan)        
    if plotmode == 'amp':
        plotvals = amp_vals
        print 'plotting amplitude solutions'
    elif plotmode == 'phase':
        plotvals = phase_vals
        print 'plotting phase solutions'
    else:
        print 'plotmode={} not defined'.format(mode)
    
    #now setup the plotting environment
    #total number of plots is number of scans
    nplots = len(scan_list)
    ny = int(np.ceil(nplots/float(nx)))
    #and set up size that I will want
    #say 3 inches per plot
    xsize = nx*plotsize
    ysize = ny*plotsize
    #and I want global limits, for best comparison
    #this could potentially change
    #will have a total offset of 
    if ymin == 0:
        ymin = np.nanmin(plotvals)
    if ymax == 0:
        ymax = np.nanmax(plotvals)
    xmin = np.nanmin(freqs)
    xmax = np.nanmax(freqs)
    
    #create figure object
    #this is what I return from program
    fig= plt.figure(figsize=(xsize,ysize))
    plt.xlim(xmin,xmax) # Limit the plot to the minimum and maximum frequencies
    plt.suptitle('Calibrated data {0}, normalization {1}'.format(plotmode,norm), fontsize='large')


    for n,(scan,beam) in enumerate(zip(scan_list,beam_list)):
        plt.subplot(ny, nx, n+1)
        for a,ant in enumerate(ant_names):
            plt.scatter(freqs[:], plotvals[a,:,pol,n],label=ant,
                       marker=',',s=0.5)
#        for sol in range(plotvals.shape[2]):
#            plt.plot(plotfreqs, (plotvals[a,:,sol] + np.full(plotvals.shape[1],offset*sol)))
        plt.title('{0}, beam {1}'.format(scan,beam))
        plt.xlim(xmin,xmax) # Limit the plot to the minimum and maximum frequencies
        plt.ylim(ymin,ymax)
           
    plt.legend()
    
    if figname != '':
        plt.savefig(figname)
    
    return fig


"""

def get_time_sequence_for_beam():
    #Get a series of scans for single beam
    
def compare_scan_solution_amp():
    #This will compare scan solutions to reference (given) for amplitude
    #will return an array that is amp solution for each scan divided by reference
    #these scans can be separated by time or beam
    
def compare_scan_solution_phase():
    #This will compare scan solutions to reference (given) for phase
    #will return an array that is phase solution for each scan divided by reference
    
def compare_scan_solution_bp_phase():
    #This will compare scan BP solutions to reference (given) for phase
    #will return a 2-D array (freq, scan) 
    #that is BP phase solution for each scan divided by reference 
        

    
def plot_compare_beam_amp():
    #this will generate plots that show comparison in amplitude solution
    #between different scans, labelled/plotted as function of beam
    
def plot_compare_beam_phase():
    #this will generate plots that show comparison in phase solution
    #between different scans, labelled/plotted as function of beam
    
def plot_compare_beam_bp_amp():
    #this will generate plots that show comparison in BP amplitude solution
    #between different scans, labelled/plotted as function of beam
    
def plot_compare_beam_bp_phase():
    #this will generate plots that show comparison in BP phase solution
    #between different scans, labelled/plotted as function of beam
    
def plot_compare_time_amp():
    #this will generate plots that show comparison in amplitude solution
    #between different scans, labelled/plotted as function of time
    
def plot_compare_time_phase():
    #this will generate plots that show comparison in phase solution
    #between different scans, labelled/plotted as function of time
    
def plot_compare_time_bp_amp():
    #this will generate plots that show comparison in BP amplitude solution
    #between different scans, labelled/plotted as function of time
    
def plot_compare_time_bp_phase():
    #this will generate plots that show comparison in BP phase solution
    #between different scans, labelled/plotted as function of time
    
"""
        