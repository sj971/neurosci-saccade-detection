# neurosci-saccade-detection
Example Matlab code for extracting useful eye position information from gaze 
coordinate data.

## About
The repository contains a Matlab script for preprocessing raw gaze coordinate 
data and detecting instances of small eye movements or 'saccades' in the time
series. The raw gaze coordinate data is loaded from an associated .mat file,
and contains gaze position data from 120 separate trials of a behavioral task.

During the task, the observer was instructed to maintain a fixed eye position; 
large deviations from stable fixation are notable on certain trials - the
positional information and velocity/acceleration profiles change abruptly in 
the traces. An example saccade trial is illustrated below, with associated 
rule-of-thumb thresholds indicated by the horizontal dashed lines.

![Saccade detection](https://raw.githubusercontent.com/sj971/neurosci-saccade-detection/sample_saccade.png “Saccade detection”)
