# EphysAnalysis_Simon

Jupyter notebooks & python scripts for LFP and MUA Analysis.

General call structure for analysis:
1. Jupyter Notebook
  --> cycleDirectory
  --> mainAnalysis
  --> analysisFunctions
  
  Additional functions in plotting, preprocessing, loadFiles, load_probe_info. 

 Written for Python 3 Version 
 
-------------------------------------------------------------------------------


MUA_analysis.py (Notebook))
    imports load_probe_info
    imports cycleDirectoryMUA.py    +3helpers
      imports mainAnalysis.py   +3helpers
        imports analysisFunctions   +3helpers

  
  We low pass filter for LFP (500Hz) extracting synchronized population activity. High pass for MUA (800), capturing rapid changes in potential, like APs. Correct? 
  We feed a uV signal into the high pass filter (which somehow filters the input uV signal based on frequency. How?) and outputs a transformed(eg no value remains at its origin amplitude?) uV signal? 
  --> The Filter selects for rapid changes (this is where the frequency comes from) and seems to output the changes (threshold-ed derivative). However, this still mismatches with a Hz (s^-1) threshold versus a uV over time (uV*s^1) signal.

  For MUA:
    8 std's as a threshold. Which distribution does the refer to? All Frequencies after highpass? 

  For LFP:
    current density analysis. Not needed for MUA?

AcuteRecording notes:

Stimulus Order: {C1, C2, D1, B1}  # just baseline to check recording location.
		{ALT 	ALTC}   # 
		{PSCC1	PSC1	O10C1	DAC1	DAC2	O10C2	O25UC1	O25C1	MS	O25C2	O25UC2	DOC1	DOC2	EQ	PSC2	PSCC2}

  How is this to be interpreted?Some paradigms seem unfamiliar.

  Flip flop not always consecutive? 

Triggers_Predeviant.dat, Triggers_Postdeviant.dat
are subsets of Triggers_Standard.dat right?

'DAC1 - deviant alone C1'
'DAC2 - deviant alone C2'
'DOC1 - deviant omission C1'
'DOC2 - deviant omission C2'
'MS   - Many Standards'
'O10C1 - Oddball 10% C1'
'O10C2 - Oddball 10% C2'
'O25C1 - Oddball 25% C1'
'O25C2 - Oddball 25% C2'
'O25UC1 - Oddball Uniform 25% C1'
'O25UC2 - Oddball Uniform 25% C2'

Analysis idea: spike filed coherence SFC

Comments in the trigger .dat file look like timepoints in seconds, however the comment says "unit=uV"?

--------------------------------------------------------------------------------


So. What data do we get from the preprocessing?

We iterate the mouse+paradigmC1/C2 (4 mice, 5C1+5C2+1MS=11)
  then, iterate the stimulus type (predev, postdev...)
    then, iterate over 32 channels. For each channel (binary uV time series file) we compute the following:
      - high-pass filtered uV in time frame of interest: -50-200ms (#stimulus, 8000 time steps)
      Each of the following for + voltage (intracell. hyperpolarization and - voltage (intracell. depolarization) ):
        - crossings: boolean array (#stimulus, -50-200ms): spikes (>8 std's)
        - timestamps: spike-timestamp array (#stimulus, max#spikesOfSomeStimulus#) - SAVED AS CSV (most dense information)
        - waveforms: time context for each spike (nested list of length = #allSpikes, inside -2ms, +2ms uV array)
        - firing rate: for 5ms bins in -50-200ms (1D 50 elements), count all spikes across stimuli  - SAVED DOWNSTREAM Triggers_Deviant_FiringRates.CSV

      On top, compute the following summary metrics:
        metadata:
        - channel filename
        - totTrials (#stimulus)
        metrics (all neg spikes)
        - #trials w/o spikes
        - Avg time to first spike
        - Sum of all spikes across trials
        - Avg spikes post stimulus 
        - slope of # spikes over trials
        - Gamma and Gaussian fits parameters, one pre stimulus, one post

      return firing rates and summary (above)
    
    2 aggregated files for 32 channels:
      1: Triggers_Deviant_FiringRates.csv
        Columns: 50 time bins, Rows: 32 channels, int counts as values

      2: Triggers_Deviant_SpikeSummary.csv:
        Columns: 17 summary stats, Rows: 32 channels, floats of different kind

---------------------------------------------------------------------------
The data:

progress:
  4/5 channels down post response at 100ms, reported by wolfger
  Deviant alone as control to subtract from.

plans:
  baseline correction
  Raster plots
  SI index
  Spike Filed Coherence 


f = feature. In my setup, relates to whisker: f1 = C1 (PW, recorded barrel), f2 = C2(AW). Imortantly, this definition is static, it doesn't change over experiments.
What does change is the assignment of a feature to 'standard' or 'deviant' , ie. the flip flip method. This should help for computing SI indices that make sense.

- implement paradigm specifc baseline subtraction
- implement paradigm specific unit/nounit 
- artifact detection through trial based correlation analysis/ look at waveforms
- Try analysis with 5 std's 
- implement variable SI index calculation 





NT Group Meeting, Wednesdays 10AM: 

https://ethz.zoom.us/j/835248730

For all further meetings between us: 

https://ethz.zoom.us/j/96761253748

SSA (Wolfger):
https://ethz.zoom.us/j/395869647

    { "key": "shift+alt+down",   "command": "editor.action.copyLinesDownAction",
        "when": "editorTextFocus && !editorReadonly" },
    { "key": "shift+alt+up",     "command": "editor.action.copyLinesUpAction",
        "when": "editorTextFocus && !editorReadonly" },


Avg_peak_neg_TS_ms? What is this?