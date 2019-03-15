# commissioning
Scripts and notebooks for commissioning

This is a common repository for working on commissiong tasks. 
The guidelines for best practices are:
- If working on existing code/notebooks, check it out and push it back when you are done actively working for the day.
- Example tutorial notebooks go directly in the commisisoning subfolder; tests run on additional multiple datasets go in the subdirectory `testing` within each commissioning folder.

The structure of this repository is: <br>
->Commissioning tasks<br>
-->python modules for the task, tutorial notebook<br>
--->testing: notebooks run on datasets to complete the commissioning tasks

## Notes on cross-calibration commissioning task
19 February 2019: Starting a massive update to rework code to rely on atdbquery and use observations in ATDB for undertaking the task. Tag the existing version as v1.0.

15 March 2019: Created functionality for looking at calibrator scans in context of 40 beam system tuning. Create plots of raw data (amplitude and phase, averaged over baselines/times for each antenna) and plots of bandpass solutions (amplitude and phase). Usage:

`python get_ccal_scans.py`

prints a list of all the sets of calibrator scans. Select a scan set, identified by first scan and source, e.g., `190314075_3C147`.

`python plot_ccal_scans.py 190314075_3C147 /data/adams/test -p --/path/to/write/images/to`
