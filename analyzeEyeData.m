function eyedat = analyzeEyeData
% function eyeDat = analyzeEyeData
%
% Output:
%   eyeDat:            data structure containing processed eye data,
%                      saccades, blink info. etc.
%
% The script runs through a block of trials (from a behavioral experiment),
% and processes trial-associated eye position traces, detecting instances
% of trials with saccades. The eye traces were originally read from an 
% Eyelink .edf file (using MGL functions, see Justin Gardner's MGL code). 
% They are provided here in .mat file format for ease of use. The
% final plotting section plots the velocity/acceleration of eye movements 
% for individual trials. Instances of saccades are obvious, with abrupt
% changes in the velocity and acceleration profiles.
%
% 
% Notes:
%
% 1. As eye movement analysis is rather time-consuming, it is a good idea
% to do the eye analysis separately from other behavioral analyses, and
% save a record of e.g., trial indices where saccades were made.
%
% 2. From the experimental design, we are interested in particular phases 
% of the trial (specifically, trial segments 4 through 7). Note the segment
% timings (in s):
%
%       task{1}.segmin = [Inf 1 .1 .6 .2 .6 .4 1.2 0.8];
%       task{1}.segmax = [Inf 1 .1 .6 .2 .6 .4 1.2 1.2];
%
% This reduced window gives ~700 position samples at 500Hz). Note that 
% initial 'inf' appears in segtimes, but with duration of only a single 
% frame or so (this is time just after observer has remained fixated for 
% 250ms and trial starts).
%
% 3. A few additional notes on analysis: (i) the median position just prior
% to stimulus onset is subtracted from subsequent position data (to limit
% possible effects of drift across session). Use of the median should limit
% the possibility of introducing artificts in shifts of position. This
% should't effect the detection of saccades/blinks or the sign of saccade
% directions; (ii) blink samples are removed towards the end of the
% pre-analysis (after velocity/acceleration are calculated but prior to
% saccade detection). This removes the possibility of blink removal edges
% creating artifacts in the velocity/acceleration calculation, while also
% ensuring that blink samples don't enter into the subsequent saccade
% detection routine, or any possible secondary analyses (e.g., of median
% position). The saccade detection routine works by searching for sequences
% in the velocity/acceleration vectors that are continuous (i.e., with no
% breaks) and above set thresholds. This should be pretty robust; (iii) the
% saccade detection routine outputs several variables for each saccade,
% including xDist and yDist (the signs of which may be useful for possible
% secondary analyses e.g., on saccade direction).
%
% 4. Finally, note that the script below is based roughly on approaches in
% papers such as Engbert & Kliegl (2003), and parts of it were adapted from
% some earlier eye analysis code of mine (Jackson et al., 2008, PLoS ONE).

%% 1. Set appropriate variables

% common parameters
nTrials = 120;
sample_window = 0.002; %2ms i.e., 500Hz tracking
nEdgeSamples = 50; %number of samples to discard either side of blink
tVelocity = 30; %threshold velocity for saccades (deg/s)
tAcceleration = 8000; %threshold acceleration for saccades (deg/s^2)
tDist = 0.5; %threshold (Euclidean) distance for saccades (deg)


%% 2. Analyze eye data

% load data block
dat = load('exampleEyeData.mat');

% analyze eye position trial by trial
for t = 1:nTrials
    
    %%%%%%% (0) %%%%%%%%%%
    
    % calculate segtimes relative to trial start
    segtimes = dat.trials(t).segtime - dat.trials(t).segtime(1);
    
    % find start and end times
    t0 = segtimes(4); %start time of first stimulus
    tEnd = segtimes(7); %end of second stimulus/start of delay
    
    % ...and sample indices closest t0/tEnd
    t0_time = abs(dat.eye.time - t0); %dat.eye.time is common to all trials
    tEnd_time = abs(dat.eye.time - tEnd);
    [~, t0_ind] = min(t0_time);
    [~, tEnd_ind] = min(tEnd_time);
    t0_ind = t0_ind - 4;
    tEnd_ind = tEnd_ind + 4;
    % ...extend window by 4 samples either side to ensure velocity
    % and acceleration calculation below covers entire stimulation
    % period (vel/acc sliding windows chop off 2 samples each
    % at either end of calculation)
    
    %%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%% (1) %%%%%%%%%%
    
    % find baseline eye position (must remove NaNs from xPos
    % samples later for calculation of median, but not for
    % subtraction from and extraction of relevant samples below)
    baselineX = dat.eye.xPos(t, 1:(t0_ind - 1));
    baselineY = dat.eye.yPos(t, 1:(t0_ind - 1));
    baselineX_NaN = isnan(baselineX);
    baselineY_NaN = isnan(baselineY);
    baselineX(baselineX_NaN) = [];
    baselineY(baselineY_NaN) = [];
    
    % get relevant eye position data aligned to baseline, and
    % also relevant pupil samples
    xPos = dat.eye.xPos(t, t0_ind:tEnd_ind) - median(baselineX);
    yPos = dat.eye.yPos(t, t0_ind:tEnd_ind) - median(baselineY);
    pupil = dat.eye.pupil(t, t0_ind:tEnd_ind);
    
    %%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%% (2) %%%%%%%%%%
    
    % blink detection/removal
    noPupil = isnan(pupil);
    eyedat.blink{t} = 0;
    
    % if trial has blink
    if max(noPupil) == 1
        
        % note instance of blink
        eyedat.blink{t} = 1;
        
        % find blink boundaries
        blinkStart = find(diff(noPupil) == 1) + 1; %(index+1 was first NaN, due to differentiation
        blinkEnd = find(diff(noPupil) == -1);
        
        % blink that started before first interval will not show up
        % in blinkStart boundary search above
        if isempty(blinkStart) == 1
            blinkStart = 1;
        end;
        % blink that ended after second interval will not show up
        % in blinkEnd boundary search above
        if isempty(blinkEnd) == 1
            blinkEnd = length(noPupil);
        end;
        
        % include samples just before/after blink,
        % accounting for timepoints that move outside analysis
        % window
        blinkStart = blinkStart - nEdgeSamples;
        blinkEnd = blinkEnd + nEdgeSamples;
        for blink = 1:length(blinkStart)
            if blinkStart(blink) < 1
                blinkStart(blink) = 1;
            end;
        end;
        for blink = 1:length(blinkEnd)
            if blinkEnd(blink) > length(noPupil)
                blinkEnd(blink) = length(noPupil);
            end;
        end;
        
        % set samples within the boundaries defined above to NaN;
        % code below should handle both single and multi-blink trials, and
        % trials in which the first or last blink extends outside
        % boundaries
        if length(blinkStart) == length(blinkEnd)
            
            for blink = 1:length(blinkStart)
                noPupil(blinkStart(blink):blinkEnd(blink)) = 1;
            end;
            
        elseif length(blinkStart) > length(blinkEnd)
            
            for blink = 1:(length(blinkStart) - 1) % all but final blink in trial
                noPupil(blinkStart(blink):blinkEnd(blink)) = 1;
            end;
            for blink = length(blinkStart) % final blink
                noPupil(blinkStart(blink):length(noPupil)) = 1;
            end;
            
        elseif length(blinkStart) < length(blinkEnd)
            
            for blink = 1 % first blink in trial
                noPupil(1:blinkEnd(blink)) = 1;
            end;
            for blink = 2:(length(blinkEnd)) % any additional blinks
                noPupil(blinkStart(blink - 1):blinkEnd(blink)) = 1;
            end;
            
        end;
        
        %...and blink samples are removed below (after velocity
        %and acceleration calculation). Removing prior to this
        %would create possibility of blink window edges being
        %accidentally marked as saccades.
        
    end;
    
    %%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%% (3) %%%%%%%%%%
    
    % calculate velocity/acceleration using 5-sample window. See
    % Engbert and Kliegl, 2003. Denominator accounts for the
    % six sample 'differences' used in numerator (i.e., n-2 to
    % n+2 = 4 samples, n-1 to n+1 = 2 samples).
    xVel = zeros(size(xPos)); yVel = zeros(size(yPos));
    for ii = 3:(size(xPos, 2) - 2) % 2 additional samples chopped off either end (see ~line 230 above)
        xVel(ii) = (xPos(ii + 2) + xPos(ii + 1) - xPos(ii - 1) - xPos(ii - 2))/(6*sample_window);
        yVel(ii) = (yPos(ii + 2) + yPos(ii + 1) - yPos(ii - 1) - yPos(ii - 2))/(6*sample_window);
    end;
    euclidVel = sqrt((xVel.*xVel) + (yVel.*yVel));
    
    xAcc = zeros(size(xPos)); yAcc = zeros(size(yPos));
    for ii = 3:(size(xVel, 2) - 2) % 2 additional samples chopped off either end (see ~line 230 above)
        xAcc(ii) = (xVel(ii + 2) + xVel(ii + 1) - xVel(ii - 1) - xVel(ii - 2))/(6*sample_window);
        yAcc(ii) = (yVel(ii + 2) + yVel(ii + 1) - yVel(ii - 1) - yVel(ii - 2))/(6*sample_window);
    end;
    euclidAcc = sqrt((xAcc.*xAcc) + (yAcc.*yAcc));
    
    % remove blink samples from pos/velocity/acceleration
    % variables, so that blink edges are not accidentally
    % marked as saccades below
    xPos(noPupil == 1) = [];
    yPos(noPupil == 1) = [];
    xVel(noPupil == 1) = [];
    yVel(noPupil == 1) = [];
    xAcc(noPupil == 1) = [];
    yAcc(noPupil == 1) = [];
    euclidVel(noPupil == 1) = [];
    euclidAcc(noPupil == 1) = [];
    
    % save to post-processed eye data structure (pupil
    % variables may now contain more samples than the others)
    eyedat.pupil{t} = pupil;
    eyedat.xPos{t} = xPos;
    eyedat.yPos{t} = yPos;
    eyedat.xVel{t} = xVel;
    eyedat.yVel{t} = yVel;
    eyedat.xAcc{t} = xAcc;
    eyedat.yAcc{t} = yAcc;
    eyedat.euclidAcc{t} = euclidAcc;
    eyedat.euclidVel{t} = euclidVel;
    eyedat.noPupil{t} = noPupil;
    
    %%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%% (4) %%%%%%%%%%
    
    % saccade detection
    eyedat.saccades{t} = [];
    candidates = find(euclidVel > tVelocity);
    if candidates
        
        % check for multiple candidate saccades in single
        % trial, using threshold parameters defined at top
        % (see Engbert & Kliegl papers, and Eyelink manual)
        saccades = [];
        diffCandidates = diff(candidates);
        breaks = [0 find(diffCandidates > 1) size(candidates, 2)];
        for jj = 1:(size(breaks, 2) - 1)
            
            % find individual candidate saccades
            saccade = [candidates(breaks(jj) + 1) candidates(breaks(jj + 1))];
            
            % exceeds acceleration threshold?
            peakAcceleration = max(euclidAcc(saccade(1):saccade(2)));
            if peakAcceleration > tAcceleration
                
                % exceeds amplitude threshold?
                xDist = xPos(saccade(2)) - xPos(saccade(1));
                yDist = yPos(saccade(2)) - yPos(saccade(1));
                euclidDist = sqrt((xDist*xDist) + (yDist*yDist));
                if euclidDist > tDist
                    
                    % store saccade info
                    peakVelocity = max(euclidVel(saccade(1):saccade(2)));
                    saccades = [saccades; saccade xDist yDist euclidDist peakVelocity];
                end
                
            end
            
        end
        
        % save saccade data
        eyedat.saccades{t} = saccades;
        
    end;
    
end;

% save to file
writeto = 'processedEyeData.mat';
save(writeto, 'eyedat')


%% 3. Trial-by-trial velocity/acceleration/position check for single block

% Saccade trials are very obvious, with the velocity and acceleration
% y-axis scales increasing dramatically beyond thresholds (30deg/s for
% velocity, 8000deg/s^2 for acceleration):
clf

for t = 1:nTrials;
    
    figure(1)
    
    subplot(3,1,1)
    plot(eyedat.euclidVel{t}, 'b', 'LineWidth', 2);
    ylabel('Velocity (deg/s)')
    ylim([0 150]);
    line([0 700], [30 30], 'LineStyle', '--', 'Color', 'r')
    
    subplot(3,1,2)
    plot(eyedat.euclidAcc{t}, 'm', 'LineWidth', 2);
    ylabel('Acceleration (deg/s^2)')
    ylim([0 15000]);
    line([0 700], [8000 8000], 'LineStyle', '--', 'Color', 'r')
    
    subplot(3,1,3)
    plot(eyedat.xPos{t}, 'k', 'LineWidth', 2);
    xlabel('Sample')
    ylabel('Position (deg)')
    ylim([-2 2]);

    pause(0.5)
    
end