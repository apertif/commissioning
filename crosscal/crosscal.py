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
sys.path.append('/home/adams/atdbquery')
from atdbquery import atdbquery
from datetime import datetime



"""Prepare data"""

def get_switching_scan_dict(maxint = 7,nskip = 1, nswitch = 30):
    #Use atdbquery to find all scans that are part of placing calibrator in (almost)every beam
    #Can run with no parameters
    #maxint is maximum integration time in minutes
    #nksip is the number of scans that can be skipped to still be part of same sequence
    #nswitch is the minimum number of total scans in a switching sequence
    allscans = atdbquery('imaging')
    scanlist=[]
    namelist=[]
    switching_scan_dict = {}
    for scan in allscans:
        starttime = scan['starttime']
        endtime = scan['endtime']
        try:
            s1 = datetime.strptime(starttime,'%Y-%m-%dT%H:%M:%SZ')
            s2 = datetime.strptime(endtime,'%Y-%m-%dT%H:%M:%SZ')
            length = s2-s1
        except TypeError:
            continue #restart loop if time is not string
        #identify short calibrator scans
        if length.seconds < maxint*60: #maxint in minutes
            scanlist.append(scan['taskID'])
            namelist.append(scan['name'])
            
    #now need to sort lists
    sorted_scan = sorted(scanlist)
    sorted_name = [name for _,name in sorted(zip(scanlist,namelist))]
    
    #then iterate through lists to parse into subsets of switching scans
    tmpscanlist=[]
    tmpnamelist=[]
    tmpobslist = []
    for scan,name in zip(sorted_scan,sorted_name):
        name_split = name.split('_') #break name down, can check if it has right structure
        #check if name has right structure, if not, skip (test obs, don't care)
        #assuming structure is Cal_Beam
        if len(name_split) == 2:
            if len(tmpscanlist) > 0: #check if already entries in list
                if ((np.abs(int(tmpscanlist[-1])-int(scan))) > (nskip+1)) or (name_split[0] != tmpnamelist[-1]):
                    #check if this scan breaks the sequence, with possiblity of nskip missed scans
                    #if source name changes, that also breaks seuqence
                    if len(tmpscanlist) >= nswitch:
                        #if so, check if long enough to add as dictionary of scans
                        keyname = tmpscanlist[0] +"_"+ tmpnamelist[0] #scan w/ calibrator
                        switching_scan_dict[keyname] = tmpobslist
                    #and reset lists if sequenc is broken
                    #scanlist is really scan plus beam
                    obs = scan+'_'+name_split[1]
                    tmpscanlist = [scan]
                    tmpobslist = [obs]
                    tmpnamelist=[name_split[0]]
                else: #if still in sequence, append
                    obs = scan+'_'+name_split[1]
                    tmpscanlist.append(scan)
                    tmpobslist.append(obs)
                    tmpnamelist.append(name_split[0])
            else: #else populate the list
                obs = scan+'_'+name_split[1]
                tmpscanlist = [scan]
                tmpobslist = [obs]
                tmpnamelist = [name_split[0]]
            
  
    #nicely format dictionary, print to screen in time order
    print 'First scan, number in switch'
    for key in sorted(switching_scan_dict):
        print '{0}  {1}'.format(key,len(switching_scan_dict[key]))
    
    return switching_scan_dict

def get_cal_data(scan_dict,basedir,scanset=None,mode='single',run=False, clearcal=False):
    #this is a wrapper function that copies, updates, flags, calibrates data
    #Use scan_dict (output of get_scan_list) to find files for copying over
    #basedir is base location to copy to
    #optional keywords:
    #scanset = None, automatically use latest scan (mode='single')
    #or all scans (mode='multiple')
    #scanset can also be a string (with mode='single') 
    #or it can be a list of strings (with mode='multiple') to plot just those scans
    #if run=False, prints data to be copied and quits

    
    #check whether or not execute commands
    if run==False:
        print 'In verification mode, will print commands to screen'
    elif run==True:
        print 'In running mode, will execute commands'
    else:
        print 'run value must be True or False'
    
    #Check the mode and scanset specification to figure out how much/what data getting
    #start with 'single' swithcing scan
    #find most recent set and get list of all - useful            
    keylist = []
    for key in sorted(scan_dict):
        keylist.append(key)
    lastset=keylist[-1]
    if mode == 'single':        
        if scanset != None:
            if scanset in keylist:
                scanset=[scanset]
            else:
                print '{} not in switching scan dict'.format(scanset)
                print 'Using last switching scan instead'
                scanset = [lastset]
        else:
            print 'Using last switching as default'
            scanset = [lastset]
    elif mode == 'multiple':
        if scanset == None:
            print 'Using all switching scans'
            scanset = keylist
        elif type(scanset) != list:
            print 'Must specify multiple scan sets in multiple mode'
            print 'Using all switching scans as default'
            scanset = keylist
        elif len(scanset) < 2:
            print 'Must specify multiple scan sets in multiple mode'
            print 'Using all switching scans as default'
            scanset = keylist
        elif all(elem in keylist for elem in scanset):
            scanset = scanset
        else:
            print 'Not all items in scanset contained in dictionary'
            print 'Using all scans by default'
            scanset = keylist
    else:
        print "mode must be 'single' or 'multiple'"
        
    #now go through data and for each scan, take appropriate action
    #record scans and beams - useful for later 
    #although I really need to change what I'm doing fundamentally here
    scanlist = []
    beamlist = []
    for scankey in scanset:
        #go through all the scans:
        for obs in scan_dict[scankey]:
            #parse scan_beam into scan,beam
            obs_split = obs.split('_')
            scan = obs_split[0]
            beam = obs_split[1]
            print 'Copying data'
            copy_scan(scan,beam,basedir,run=run)
            scanlist.append(scan)
            beamlist.append(beam)
            if run==False:
                print 'In verification mode, no data copied, will skip remaining steps'
            else:
                print 'Fixing source names'
                fix_source_name(scan,beam,basedir)
                print 'Flagging data'
                flag_scan(scan,beam,basedir)
                print 'Calibrating data'
                calibrate_scan(scan,beam,basedir)
                
    return scanlist,beamlist
                
            
    

def copy_scan(scan,beam,basedir,run=False):
    #Copy a scan/beam combination over
   
    #move into targetdir
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
    string_command = "python /home/adams/altadata/getdata_alta.py {0} {1}-{1} {2:0>2}-{2:0>2}".format(scandate,scannum,beam)
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

                
def fix_source_name(scan,beam,basedir):
    #Source names have beam in them (critical for get_scans!)
    #But that will cause calibration to crash
    #CASA isn't very smart
    #So I'll have to be smart instead and update everything
    targetdir = '{0}/{1}/00/raw'.format(basedir,scan)
    msfile = "{0}/WSRTA{1}_B{2:0>3}.MS".format(targetdir,scan,beam)
    t_field = pt.table(msfile+"::FIELD", readonly=False)
    name = t_field[0]['NAME']
    name_split = name.split('_') #split by underscore - works also if no underscore
    t_field.putcell("NAME", 0, name_split[0])  #update source name to anything before first underscore (or original name if no underscore)
    t_field.flush()
    
    
def flag_scan(scan,beam,basedir,edges=True,ghosts=True):
    preflag = apercal.preflag()
    preflag.preflag_edges = edges #option to turn on/off
    preflag.preflag_ghosts = ghosts #option to turn on/off
    #set source and pol cal to False - will iterate through updating fluxcal and running preflag one at a time
    preflag.preflag_manualflag_polcal = False
    preflag.preflag_manualflag_target = False
    preflag.preflag_aoflagger_polcal = False 
    preflag.preflag_aoflagger_target = False
    #should only find central beam, but specify anyway
    preflag.preflag_aoflagger_targetbeams = '00' 
    preflag.fluxcal = "WSRTA{0}_B{1:0>3}.MS".format(scan,beam)
    preflag.basedir = "{0}/{1}/".format(basedir,scan) #basedir for Apercal is different than my basedir, needs trailing slash
    print "Setting fluxcal to WSRTA{0}_B{1:0>3}.MS".format(scan,beam)
    #and run preflag
    print "Flagging data set {2}00/raw/WSRTA{0}_B{1:0>3}.MS".format(scan,beam,preflag.basedir)
    preflag.go()

    
def calibrate_scan(scan,beam,basedir):
    #load ccal module
    ccal = apercal.ccal()
    ccal.crosscal_transfer_to_target = False #there is no target
    ccal.fluxcal = "WSRTA{0}_B{1:0>3}.MS".format(scan,beam)
    ccal.basedir = "{0}/{1}/".format(basedir,scan) #basedir for Apercal is different than my basedir, plus trailing slash
    print "Setting fluxcal to WSRTA{0}_B{1:0>3}.MS".format(scan,beam)
    print "Calibrating data set {2}00/raw/WSRTA{0}_B{1:0>3}.MS".format(scan,beam,ccal.basedir)
    ccal.go()



"""Object specification"""    
"""
Define object classes for holding data related to scans
"""

class ScanData(object):
    #Initilailze with source name, scalist and beamlist
    #and place holders for phase and amplitude
    def __init__(self,source,basedir,scanlist,beamlist):
        self.source = source
        self.scanlist = scanlist
        self.beamlist = beamlist
        self.basedir=basedir
        self.phase = np.empty(len(scanlist),dtype=np.ndarray)
        self.amp = np.empty(len(scanlist),dtype=np.ndarray)
        
class BPSols(ScanData):
    def __init__(self,source,basedir,scanlist,beamlist):
        ScanData.__init__(self,source,basedir,scanlist,beamlist)
        self.ants = np.empty(len(scanlist),dtype=np.object)
        self.time = np.empty(len(scanlist),dtype=np.ndarray)
        self.freq = np.empty(len(scanlist),dtype=np.ndarray)
        self.flags = np.empty(len(scanlist),dtype=np.ndarray)
        self.amps_norm = np.empty(len(scanlist))
        self.phases_norm = np.empty(len(scanlist))
    
    def get_data(self):
        #get the data
        for i, (scan,beam) in enumerate(zip(self.scanlist,self.beamlist)):
            bptable = "{0}/{1}/00/raw/WSRTA{1}_B{2:0>3}.Bscan".format(self.basedir,scan,beam)
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
            
            #check for flags and mask
            amp_sols[flags] = np.nan
            phase_sols[flags] = np.nan
            
            self.ants[i] = ant_names
            self.time[i] = times
            self.phase[i] = phase_sols *180./np.pi #put into degrees
            self.amp[i] = amp_sols
            self.flags[i] = flags
            self.freq[i] = freqs
            
    def plot_amp(self):
        #plot amplitude, one plot per antenna
        #put plots in default place w/ default name
        ant_names = self.ants[0]
        #figlist = ['fig_'+str(i) for i in range(len(ant_names))]
        for a,ant in enumerate(ant_names):
            #iterate through antennas
            #set up for 8x5 plots (40 beams)
            nx = 8
            ny = 5
            xsize = nx*4
            ysize = ny*4
            plt.figure(figsize=(xsize,ysize))
            plt.suptitle('Bandpass amplitude for Antenna {0}'.format(ant))
            
            for n,(scan,beam) in enumerate(zip(self.scanlist,self.beamlist)):
                beamnum = int(beam)
                plt.subplot(ny, nx, beamnum+1)
                plt.scatter(self.freq[n][0,:],self.amp[n][a,:,0],
                            label='XX, {0}'.format(self.time[n][a]),
                            marker=',',s=1)
                plt.scatter(self.freq[n][0,:],self.amp[n][a,:,1],
                            label='YY, {0}'.format(self.time[n][a]),
                            marker=',',s=1)
                plt.title('Beam {0}'.format(beam))
                plt.ylim(0,1.8)
            plt.legend()
            plt.savefig('/home/adams/commissioning/crosscal/img/BP_amp_{0}_{1}.png'.format(ant,self.scanlist[0][0:6]))
            
    def plot_phase(self):
        #plot phase, one plot per antenna
        ant_names = self.ants[0]
        #figlist = ['fig_'+str(i) for i in range(len(ant_names))]
        for a,ant in enumerate(ant_names):
            #iterate through antennas
            #set up for 8x5 plots (40 beams)
            nx = 8
            ny = 5
            xsize = nx*4
            ysize = ny*4
            plt.figure(figsize=(xsize,ysize))
            plt.suptitle('Bandpass phases for Antenna {0}'.format(ant))
            
            for n,(scan,beam) in enumerate(zip(self.scanlist,self.beamlist)):
                beamnum = int(beam)
                plt.subplot(ny, nx, beamnum+1)
                plt.scatter(self.freq[n][0,:],self.phase[n][a,:,0],
                            label='XX, {0}'.format(self.time[n][a]),
                            marker=',',s=1)
                plt.scatter(self.freq[n][0,:],self.phase[n][a,:,1],
                            label='YY, {0}'.format(self.time[n][a]),
                            marker=',',s=1)
                plt.title('Beam {0}'.format(beam))
                plt.ylim(-180,180)
            plt.legend()
            plt.savefig('/home/adams/commissioning/crosscal/img/BP_phase_{0}_{1}.png'.format(ant,self.scanlist[0][0:6]))
            
            
        
class GainSols(ScanData):
    def __init__(self,source,basedir,scanlist,beamlist):
        ScanData.__init__(self,source,basedir,scanlist,beamlist)
        self.ants = np.empty(len(scanlist),dtype=np.object)
        self.time = np.empty(len(scanlist),dtype=np.ndarray)
        self.flags = np.empty(len(scanlist),dtype=np.ndarray)
        self.amps_norm = np.empty(len(scanlist),dtype=np.ndarray)
        self.phases_norm = np.empty(len(scanlist),dtype=np.ndarray)
        
    def get_data(self):
        for i, (scan,beam) in enumerate(zip(self.scanlist,self.beamlist)):
            gaintable = "{0}/{1}/00/raw/WSRTA{1}_B{2:0>3}.G1ap".format(self.basedir,scan,beam)
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
                
            #check for flags and mask
            amp_ant_array[flags_ant_array] = np.nan
            phase_ant_array[flags_ant_array] = np.nan
            
            self.amp[i] = amp_ant_array
            self.phase[i] = phase_ant_array * 180./np.pi #put into degrees
            self.ants[i] = ant_names
            self.time[i] = times
            self.flags[i] = flags_ant_array
            
    def plot_amp(self):
        #plot amplitude, one plot per antenna
        #put plots in default place w/ default name
        ant_names = self.ants[0]
        #figlist = ['fig_'+str(i) for i in range(len(ant_names))]
        for a,ant in enumerate(ant_names):
            #iterate through antennas
            #set up for 8x5 plots (40 beams)
            nx = 8
            ny = 5
            xsize = nx*4
            ysize = ny*4
            plt.figure(figsize=(xsize,ysize))
            plt.suptitle('Gain amplitude for Antenna {0}'.format(ant))
            
            for n,(scan,beam) in enumerate(zip(self.scanlist,self.beamlist)):
                beamnum = int(beam)
                plt.subplot(ny, nx, beamnum+1)
                plt.plot(self.time[n],self.amp[n][a,:,0],
                         label='XX, {0}'.format(self.time[n][0]))
                plt.plot(self.time[n],self.amp[n][a,:,1],
                         label='YY, {0}'.format(self.time[n][0]))
                plt.title('Beam {0}'.format(beam))
                plt.ylim(10,30)
            plt.legend()
            plt.savefig('/home/adams/commissioning/crosscal/img/Gain_amp_{0}_{1}.png'.format(ant,self.scanlist[0][0:6]))

        
class ModelData(ScanData):
    def __init__(self,source,basedir,scanlist,beamlist):
        ScanData.__init__(self,source,basedir,scanlist,beamlist)
        self.freq = np.empty(len(scanlist))
        
class CorrectedData(ScanData):
    def __init__(self,source,basedir,scanlist,beamlist):
        ScanData.__init__(self,source,basedir,scanlist,beamlist)
        self.ants = np.empty(len(scanlist))
        self.freq = np.empty(len(scanlist))

    
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

def compare_scan_solution_bp(scan_list,beam_list,basedir,norm=True,refscan=''):
    #This will compare scan BP solutions to reference (given) for amp & phase 
    #(easiest to do together at once)
    #will return a multi-d array (ant, freq, pol, scan)
    #If norm=True, amp solution is divided by reference and phase ref is subtracted from sol
    #otherwise it returns all the solutions into array
    
    
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


def plot_compare_bp_beam(scan_list,beam_list,basedir,norm=True,
                         refscan='',plotmode='amp',pol=0,nx=3,ymin=0,ymax=0,plotsize=4,
                        figname=''):
    #this will generate plots that compare BP solutions between beams
    #It can run with option amplitude or phase, default amp
    #Defaults to showing pol 0, can also change
    
   
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
        
    print np.nanmin(bp_amp_vals),np.nanmax(bp_amp_vals)
    print np.nanmin(bp_phase_vals),np.nanmax(bp_phase_vals)
    
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
    
    if figname == '':
        return fig
    else:
        plt.savefig(figname)
        


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


def compare_scan_solution_gain(scan_list,beam_list,basedir,norm=True,refscan=''):
    #This will collect all gain solutions for a scan of beams
    #If set, will nromalize to a reference.
    #Note that each scan will have different times, 
    #so will compute a scalar average
    #and use that for comparison
    #for now, assume separated by beam but by time should also work
    
    
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


def plot_compare_gain_beam(scan_list,beam_list,basedir,
                           norm=True,refscan='',plotmode='amp',
                           pol=0,nx=8,ymin=0,ymax=0,plotsize=4,
                          figname=''):
    #this will generate plots that compare BP solutions between beams
    #It can run with option amplitude or phase, default amp
    #Defaults to showing pol 0, can also change
    
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
    
    if figname == '':
        return fig
    else:
        plt.savefig(figname)


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

def compare_scan_model(scan_list,beam_list,basedir):
    #take a scan object and get the model for each observation
    #this will confirm that sources are properly named
    
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

def plot_compare_scan_model(scan_list,beam_list,basedir,plotmode='amp',
                            pol=0,nx=8,ymin=0,ymax=0,plotsize=4,
                           figname=''):
   
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
        
    if figname == '':
        return fig
    else:
        plt.savefig(figname)
        

    
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


def compare_scan_calibrated_data(scan_list,beam_list,basedir,norm=True,refscan=''):
    #can normalize to reference
    #useful to see if there are systematic issues, but wouldn't expect taht
    
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
            

def plot_compare_calibrated_data_beam(scan_list,beam_list,basedir,norm=True,
                                      refscan='',plotmode='amp',pol=0,nx=3,ymin=0,ymax=0,plotsize=4,
                                     figname=''):
    #this will generate plots that compare BP solutions between beams
    #It can run with option amplitude or phase, default amp
    #Defaults to showing pol 0, can also change

    ant_names,freqs,amp_vals,phase_vals = compare_scan_calibrated_data(scans,obsrecordfile,basedir,
                                                                       norm=norm,refscan=refscan)        
    print np.min(amp_vals),np.max(amp_vals)
    
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
    
    if figname == '':
        return fig
    else:
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
        