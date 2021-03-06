# neurosci-saccade-detection
Example Matlab code for extracting useful eye position information from gaze 
coordinate data.

## About
The repository contains a Matlab script for preprocessing raw gaze coordinate 
data and detecting instances of small eye movements or **saccades** in the time
series. This procedure is often carried out in order to remove from analysis trials 
in which an observer was not attending to a given experimental task; however, the same
routine can be modified to study how different task or stimulus parameters might affect
gaze position or eye movement dynamics.

The raw gaze coordinate data (500Hz) is loaded from an associated .mat file, and 
contains gaze position data from 120 separate trials of a behavioral task. A number
of **initial preprocessing steps** are first taken e.g., application of a median sliding 
window across the traces to remove very high frequency fluctuations, removal of missing 
data or 'blink' samples. Instances of saccade trials are then detected by applying 
rule-of-thumb velocity and acceleration thresholds to the traces. 

During the experimental task, the observer was instructed to maintain a fixed eye position; 
yet, **large deviations from stable fixation** are notable on certain trials - the gaze position 
and velocity/acceleration profiles change abruptly. An example saccade trial is illustrated 
below (threshold cut-offs are indicated by dashed red lines).

### Further details

Details of the behavioral experiments and control eyetracking data can be found in Chapter 2 of my PhD thesis here:                              
http://sj971.github.io/docs/thesis_sjackson.pdf

![Saccade detection](sample_saccade.png)
