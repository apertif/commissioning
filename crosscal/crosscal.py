#python "module" for tests and inventories related to
#confirming the ability and calculating appropriate scaling
#for transfering solutions from central beam to other beams

#load necessary packages
import os
import numpy as np
from astropy.io import ascii
import sys
sys.path.append('/home/adams/altadata')

#define a class that is used for specifying the scans to be investigated

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

            



def get_scan_list(scans,obsrecordfile):
    #This takes a ScanObject that specifies the scans wanted
    #Use this information to produce a scan and beam list
    #All data is on ALTA now, so need to search there
    #ATDB is not yet up to the task, so use Apertif Observation Record
    #provide this as a specific input in csv format
    
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
        scan_list = np.empty(int(scans.nscan),dtype=object) #set scanlist to be lenght of nscan
        beam_list = np.empty(int(scans.nscan),dtype=object) #set beamlist to length of nscan
        #iterate through each scan and add it to the scan_list
        #parse obsrecord to find relevant beam aand add beam to beamlist
        for n in xrange(int(scans.nscan)):
            scannumber = int(scans.startscan)+n #get scan number as an int
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
    
    
    return mode,scan_list,beam_list
 
  
def copy_scans(scans,obsrecordfile,targetdir,run=False):
    #This will use the output of get_scan_list 
    #to retrieve data from ALTA and move to targetdest
    #I will take same input as get_scan_list and run it here
    #Want to minimize number of calls by user
    #Maybe the proper way to do this is to have global variables, 
    #But that's a step too complicated for me, I think
    mode,scan_list,beam_list = get_scan_list(scans,obsrecordfile)
    #get the values I need
    
    #get the data
    
    #first, move into targetdest (create if needed)
    #assumption is all data goes in same directory
    if os.path.exists(targetdir):
        pass
    else:
        os.makedirs(targetdir)
    
    #then change to working directory to copy data
    #think this is how the usage of ALTA data module from Vanessa is
    os.chdir(targetdir)
    print 'Moved to and copying data to {0}'.format(os.getcwd())
    
    if run==False:
        print 'In verification mode, will print commands to screen'
    elif run==True:
        print 'In running mode, will execute commands'
    else:
        print 'Mode is not defined. Choose verify or run'
    
    #now iterate through scan/beam lists and run command to retrieve data
    for scan,beam in zip(scan_list,beam_list):
        #first parse scan as required:
        scandate = scan[0:6]
        scannum = scan[6:9]
        #format string command to match usage:
            # ALTA data transfer: Uses the iROD client to transfer data from ALTA
            # Example usage: >> python getdata_alta.py 180316 004-010 00-36
            # V.A. Moss (vmoss.astro@gmail.com)
        #hope is that since I did sys.path.append('/home/adams/altadata'),
        #I don't have to specify full path
        string_command = "python getdata_alta.py {0} {1} {2:0>2}".format(scandate,scannum,beam)
        #!!!!!!!Important question - do I have specify beams with double integers?
        #Thayt's annoying because it doesn't match naming scheme
        #But maybe I can do with format - yep!
        if run==False:
            print string_command
        elif run==True:
            #and run the command:
            pass #just until check with Vanessa
            #os.system(string_command)
    
    
"""
def flag_scans():
    #This will use output of get_scan_list
    #plus apercal.preflag() to flag the data
    #Note that this likely will not work until the preflag error is fixed
    #So I probably won't use it at first    

def convert_scans():
    #Convert to miriad format
    
def calibrate_scans():
    #This iwll use scan_list,beam_list from get_scan_list
    #plus apercal.ccal to calibrate the data
    #HAVE TO CHANGE SOURCE NAME!
    


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
        
def compare_scan_solution_bp_amp():
    #This will compare scan BP solutions to reference (given) for amp
    #will return a 2-D array (freq, scan)
    #that is BP amp solution for each scan divided by reference
    
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
        