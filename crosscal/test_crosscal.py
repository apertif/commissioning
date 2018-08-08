#collection of tests for crosscal.py

#I will create a separate TestObject for each function in crosscal.py
#That way, I can just run tests on a single object if I only update one function

import sys
sys.path.append('/home/adams/apertif/CROSSCAL')
import crosscal as cc

class TestClass_ScanSpecification(object):
    #define scanspecification as a class
    #with ability to set different attributes
    #all initialized to empty strings
    def test_initialize_startdate_ScanSpec(self):
        scans = cc.ScanSpecification()
        assert scans.startdate == ''
    def test_initialize_enddate_ScanSpec(self):
        scans = cc.ScanSpecification()
        assert scans.enddate == ''
    def test_initialize_beam_ScanSpec(self):
        scans = cc.ScanSpecification()
        assert scans.beam == ''        
    def test_initialize_startscan_ScanSpec(self):
        scans = cc.ScanSpecification()
        assert scans.startscan == ''
    def test_initialize_endscan_ScanSpec(self):
        scans = cc.ScanSpecification()
        assert scans.endscan == ''        
    def test_initialize_nscan_ScanSpec(self):
        scans = cc.ScanSpecification()
        assert scans.nscan == ''  
    def test_set_bad_string_startdate(self):
        scans=cc.ScanSpecification()
        scans.setstartdate('nd39sl2k')
        assert scans.startdate == ''
    def test_set_short_startdate(self):
        scans =cc.ScanSpecification()
        scans.setstartdate('180512')
        assert scans.startdate == ''
    def test_set_early_startdate(self):
        scans=cc.ScanSpecification()
        scans.setstartdate('20050812')
        assert scans.startdate == ''
    def test_set_good_startdate(self):
        scans=cc.ScanSpecification()
        scans.setstartdate('20180512')
        assert  scans.startdate ==  '20180512'
    def test_set_long_string_startdate(self):
        scans =cc.ScanSpecification()
        scans.setstartdate('2018asdkglbnasldkgjaesltkj')
        assert scans.startdate == ''
    def test_set_number_startdate(self):
        scans =cc.ScanSpecification()
        scans.setstartdate(20180512)
        assert scans.startdate == '20180512'
    def test_set_bad_string_enddate(self):
        scans=cc.ScanSpecification()
        scans.setenddate('nd39sl2k')
        assert scans.enddate == ''
    def test_set_short_enddate(self):
        scans =cc.ScanSpecification()
        scans.setenddate('180512')
        assert scans.enddate == ''
    def test_set_early_enddate(self):
        scans=cc.ScanSpecification()
        scans.setenddate('20050812')
        assert scans.enddate == ''
    def test_set_good_enddate(self):
        scans=cc.ScanSpecification()
        scans.setenddate('20180512')
        assert  scans.enddate ==  '20180512'
    def test_set_long_string_enddate(self):
        scans =cc.ScanSpecification()
        scans.setenddate('2018asdkglbnasldkgjaesltkj')
        assert scans.enddate == ''        
    def test_set_number_wenddate(self):
        scans =cc.ScanSpecification()
        scans.setenddate(20180512)
        assert scans.enddate == '20180512'
    def test_set_beam_good_string(self):
        scans = cc.ScanSpecification()
        scans.setbeam('03')
        assert scans.beam == '03'
    def test_set_beam_single_int_string(self):
        scans = cc.ScanSpecification()
        scans.setbeam('3')
        assert scans.beam == '03'
    def test_set_beam_single_integer(self):
        scans = cc.ScanSpecification()
        scans.setbeam(3)
        assert scans.beam == '03'
    def test_set_beam_double_integer(self):
        scans = cc.ScanSpecification()
        scans.setbeam(15)
        assert scans.beam == '15'
    def test_set_beam_too_large_string(self):
        scans = cc.ScanSpecification()
        scans.setbeam('42')
        assert scans.beam == ''
    def test_set_beam_too_large_int(self):
        scans = cc.ScanSpecification()
        scans.setbeam(42)
        assert scans.beam == ''
    def test_set_beam_too_small_int(self):
        scans = cc.ScanSpecification()
        scans.setbeam(-1)
        assert scans.beam == ''
    def test_set_beam_too_small_string(self):
        scans = cc.ScanSpecification()
        scans.setbeam(-1)
        assert scans.beam == ''
    def test_set_beam_float_string(self):
        scans = cc.ScanSpecification()
        scans.setbeam('21.5')
        assert scans.beam == ''
    def test_set_short_startscan(self):
        scans =cc.ScanSpecification()
        scans.setstartscan('18051201')
        assert scans.startscan == ''
    def test_set_early_startscan(self):
        scans=cc.ScanSpecification()
        scans.setstartscan('080512001')
        assert scans.startscan == ''
    def test_set_good_startscan(self):
        scans=cc.ScanSpecification()
        scans.setstartscan('180512001')
        assert  scans.startscan ==  '180512001'
    def test_set_long_string_startscan(self):
        scans =cc.ScanSpecification()
        scans.setstartscan('18asdkglbnasldkgjaesltkj')
        assert scans.startscan == ''
    def test_set_number_startscan(self):
        scans =cc.ScanSpecification()
        scans.setstartscan(180512001)
        assert scans.startscan == '180512001'        
    def test_set_short_endscan(self):
        scans =cc.ScanSpecification()
        scans.setendscan('18051201')
        assert scans.endscan == ''
    def test_set_early_endscan(self):
        scans=cc.ScanSpecification()
        scans.setendscan('080512001')
        assert scans.endscan == ''
    def test_set_good_endscan(self):
        scans=cc.ScanSpecification()
        scans.setendscan('180512001')
        assert  scans.endscan ==  '180512001'
    def test_set_long_string_endscan(self):
        scans =cc.ScanSpecification()
        scans.setendscan('18asdkglbnasldkgjaesltkj')
        assert scans.endscan == ''
    def test_set_number_endscan(self):
        scans =cc.ScanSpecification()
        scans.setendscan(180512001)
        assert scans.endscan == '180512001'         
    def test_set_nscan_good_int(self):
        scans = cc.ScanSpecification()
        scans.setnscan(235)
        assert scans.nscan == '235'
    def test_set_nscan_good_string(self):
        scans = cc.ScanSpecification()
        scans.setnscan('235')
        assert scans.nscan == '235'
    def test_set_nscan_float(self):
        scans=cc.ScanSpecification()
        scans.setnscan(235.7)
        assert scans.nscan == '235'
    def test_set_nscan_float_string(self):
        scans = cc.ScanSpecification()
        scans.setnscan('235.7')
        assert scans.nscan == '235'
    def test_set_nscan_bad_string(self):
        scans = cc.ScanSpecification()
        scans.setnscan('asgase')
        assert scans.nscan == ''
  
        
class TestClass_get_scan_list(object):
    def test_get_scan_list_startdateonly(self):
        scans = cc.ScanSpecification()
        scans.setstartdate('20180504')
        mode,scanlist,beamlist = cc.get_scan_list(scans)
        assert mode == None
    def test_get_scanlist_enddateonly(self):
        scans = cc.ScanSpecification()
        scans.setenddate('20180507')
        mode,scanlist,beamlist=cc.get_scan_list(scans)
        assert mode == None
    def test_get_scanlist_startdate_enddate_nobeam(self):
        scans=cc.ScanSpecification()
        scans.setstartdate('20180403')
        scans.setenddate('20180508')
        mode,scanlist,beamlist=cc.get_scan_list(scans)
        assert mode == None
    def test_get_scanlist_beamonly(self):
        scans =cc.ScanSpecification()
        scans.setbeam('23')
        mode,scanlist,beamlist=cc.get_scan_list(scans)
        assert mode == None
    def test_get_scanlist_variabilitymode(self):
        scans = cc.ScanSpecification()
        scans.setbeam('23')
        scans.setenddate('20180607')
        scans.setstartdate('20180504')
        mode,scanlist,beamlist=cc.get_scan_list(scans)
        assert mode == 'variability'
    def test_get_scanlist_startscanonly(self):
        scans = cc.ScanSpecification()
        scans.setstartscan('180403003')
        mode,scanlist,beamlist = cc.get_scan_list(scans)
        assert mode == None
    def test_get_scanlist_startendscan(self):
        scans = cc.ScanSpecification()
        scans.setstartscan('180403001')
        scans.setendscan('180403032')
        mode,scanlist,beamlist = cc.get_scan_list(scans)
        assert mode == 'switch'
    def test_get_scanlist_startscan_nscan(self):
        scans = cc.ScanSpecification()
        scans.setstartscan('180403001')
        scans.setnscan(23)
        mode,scanlist,beamlist = cc.get_scan_list(scans)
        assert mode == 'switch'
    def test_get_scanlist_nscanonly(self):
        scans = cc.ScanSpecification()
        scans.setnscan(24)
        mode,scanlist,beamlist=cc.get_scan_list(scans)
        assert mode == None
    def test_get_scanlist_endscanonly(self):
        scans = cc.ScanSpecification()
        scans.setendscan('180403023')
        mode,scanlist,beamlist = cc.get_scan_list(scans)
        assert mode == None
        