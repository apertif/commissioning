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
            
    def plot_phase(self):
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
            plt.suptitle('Gain phase for Antenna {0}'.format(ant))
            
            for n,(scan,beam) in enumerate(zip(self.scanlist,self.beamlist)):
                beamnum = int(beam)
                plt.subplot(ny, nx, beamnum+1)
                plt.plot(self.time[n],self.phase[n][a,:,0],
                         label='XX, {0}'.format(self.time[n][0]))
                plt.plot(self.time[n],self.phase[n][a,:,1],
                         label='YY, {0}'.format(self.time[n][0]))
                plt.title('Beam {0}'.format(beam))
                plt.ylim(10,30)
            plt.legend()
            plt.savefig('/home/adams/commissioning/crosscal/img/Gain_phase_{0}_{1}.png'.format(ant,self.scanlist[0][0:6]))

        
class ModelData(ScanData):
    def __init__(self,source,basedir,scanlist,beamlist):
        ScanData.__init__(self,source,basedir,scanlist,beamlist)
        self.freq = np.empty(len(scanlist),dytpe=np.ndarray)
        
    def get_data(self):
        for i, (scan,beam) in enumerate(zip(self.scanlist,self.beamlist)):
            msfile = "{0}/{1}/00/raw/WSRTA{1}_B{2:0>3}.MS".format(self.basedir,scan,beam)
            taql_command = "SELECT abs(gmeans(MODEL_DATA)) AS amp, arg(gmeans(MODEL_DATA)) AS phase FROM {0}".format(msfile)
            t = pt.taql(taql_command)
            amp = t.getcol('amp')[0,:,:]
            phase = t.getcol('phase')[0,:,:]
            taql_freq = "SELECT CHAN_FREQ FROM {0}::SPECTRAL_WINDOW".format(msfile)
            t = pt.taql(taql_freq)
            freqs = t.getcol('CHAN_FREQ')[0,:]
            
            self.amp[i] = amp
            self.phase[i] = phase
            self.freq[i] = freqs
            
    def plot_amp(self):
        #plot amplitude, one subplot per beam
        #put plots in default place w/ default name
        nx = 8
        ny = 5
        xsize = nx*4
        ysize = ny*4
        plt.figure(figsize=(xsize,ysize))           
        plt.suptitle('Model amplitude')
            
        for n,(scan,beam) in enumerate(zip(self.scanlist,self.beamlist)):
            beamnum = int(beam)
            plt.subplot(ny, nx, beamnum+1)
            plt.plot(self.freq[n],self.amp[n][:,0],
                     label='XX')
            plt.plot(self.freq[n],self.amp[n][:,1],
                     label='YY')
            plt.title('Beam {0}'.format(beam))
            #plt.ylim(10,30)
            plt.legend()
            plt.savefig('/home/adams/commissioning/crosscal/img/Model_amp_{0}.png'.format(self.scanlist[0][0:6]))
            
    def plot_phase(self):
        #plot amplitude, one subplot per beam
        #put plots in default place w/ default name
        nx = 8
        ny = 5
        xsize = nx*4
        ysize = ny*4
        plt.figure(figsize=(xsize,ysize))           
        plt.suptitle('Model phase')
            
        for n,(scan,beam) in enumerate(zip(self.scanlist,self.beamlist)):
            beamnum = int(beam)
            plt.subplot(ny, nx, beamnum+1)
            plt.plot(self.freq[n],self.phase[n][:,0],
                     label='XX')
            plt.plot(self.freq[n],self.phase[n][:,1],
                     label='YY')
            plt.title('Beam {0}'.format(beam))
            #plt.ylim(10,30)
            plt.legend()
            plt.savefig('/home/adams/commissioning/crosscal/img/Model_phase_{0}.png'.format(self.scanlist[0][0:6]))
            
        
class CorrectedData(ScanData):
    def __init__(self,source,basedir,scanlist,beamlist):
        ScanData.__init__(self,source,basedir,scanlist,beamlist)
        self.ants = np.empty(len(scanlist),dytpe=np.object)
        self.freq = np.empty(len(scanlist),dtype=np.array)
        
    def get_data(self):
        for i, (scan,beam) in enumerate(zip(self.scanlist,self.beamlist)):
            msfile = "{0}/{1}/00/raw/WSRTA{1}_B{2:0>3}.MS".format(self.basedir,scan,beam)
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
                                "arg(gmeans(CORRECTED_DATA[FLAG])) AS phase FROM {0} "
                                "WHERE ANTENNA1!=ANTENNA2 && "
                                "(ANTENNA1={1} || ANTENNA2={1})").format(msfile,ant)
                t = pt.taql(taql_command)
                test=t.getcol('amp')
                amp_ant_array[ant,:,:] = t.getcol('amp')[0,:,:]
                phase_ant_array[ant,:,:] = t.getcol('phase')[0,:,:]
                
            self.phase[i] = phase_ant_array
            self.amp[i] = amp_ant_array
            self.freq[i] = freqs
            self.ants[i] = ant_names
            
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
            plt.suptitle('Corrected amplitude for Antenna {0}'.format(ant))
            
            for n,(scan,beam) in enumerate(zip(self.scanlist,self.beamlist)):
                beamnum = int(beam)
                plt.subplot(ny, nx, beamnum+1)
                plt.plot(self.freq[n],self.amp[n][a,:,0],
                         label='XX, {0}'.format(scan))
                plt.plot(self.freq[n],self.amp[n][a,:,3],
                         label='YY, {0}'.format(scan))
                plt.title('Beam {0}'.format(beam))
                plt.ylim(10,30)
            plt.legend()
            plt.savefig('/home/adams/commissioning/crosscal/img/Corrected_amp_{0}_{1}.png'.format(ant,self.scanlist[0][0:6]))
            
    def plot_phase(self):
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
            plt.suptitle('Corrected phase for Antenna {0}'.format(ant))
            
            for n,(scan,beam) in enumerate(zip(self.scanlist,self.beamlist)):
                beamnum = int(beam)
                plt.subplot(ny, nx, beamnum+1)
                plt.plot(self.freq[n],self.phase[n][a,:,0],
                         label='XX, {0}'.format(scan))
                plt.plot(self.freq[n],self.phase[n][a,:,3],
                         label='YY, {0}'.format(scan))
                plt.title('Beam {0}'.format(beam))
                plt.ylim(10,30)
            plt.legend()
            plt.savefig('/home/adams/commissioning/crosscal/img/Corrected_phase_{0}_{1}.png'.format(ant,self.scanlist[0][0:6]))
                         

    



