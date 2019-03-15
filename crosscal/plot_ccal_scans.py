"""
Script to automatically run crosscal plots
Requires a scan number
Optionally takes a directory for writing plots
"""

#import os
import crosscal as cc
import argparse
from timeit import default_timer as timer 



parser = argparse.ArgumentParser(description='Generate plots of calibrator scans')

# 1st argument: File name
parser.add_argument("scanset", help='identifier for a scanset, found via get_ccal_scans')
parser.add_argument("basedir",help='Location to put data')
#parser.add_argument("fluxcal", help='Fluxcal name')

parser.add_argument('-p', '--path',default=None,
                    help='Destination for images') 

args = parser.parse_args()


split_scanset = args.scanset.split('_')

cal = split_scanset[1]
scan = split_scanset[0]

start = timer()
#Get data
print 'Getting data and doing flagging/calibration'
scandict = cc.get_switching_scan_dict()
scanlist,beamlist = cc.get_cal_data(scandict,args.basedir,run=True,mode='single',scanset=args.scanset)
end = timer()

print 'Elapsed time to copy, flag and calibrate data is {} minutes'.format((end - start)/60.) 


#get BP solution plots
start = timer()
print 'getting bandpass solution plots'
BP = cc.BPSols(cal,args.basedir,scanlist,beamlist)
BP.get_data()
BP.plot_amp(imagepath=args.path)
BP.plot_phase(imagepath=args.path)

print 'Done with bandpass plots'


print 'getting raw data'

#Get Raw data
Raw = ccplots.RawData(cal,args.basedir,scanlist,beamlist)
Raw.get_data()
Raw.plot_amp(imagepath=args.path)
Raw.plot_phase(imagepath=args.path)

print 'Done with plotting raw data'


end = timer()
print 'Elapsed time to generate calibration scan plots is {} minutes'.format((end - start)/60.) 
#time in minutes


