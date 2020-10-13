function [baseline, doPupil] = get_baseline(b, t, doPupil)
% b = block generated by generateGenBlock
% t as loaded by loadSVDfiles; so timepoints from imaging
% or from load_ROI_lfp for timepoints from ephys

if nargin < 3
    doPupil = false;
end

SetDefaultDirs;
% BehDir   = DIRS.expInfo;
% BehDir = '\\zclone.cortexlab.net\Data\expInfo';
BehDir = '\\zserver.cortexlab.net\Data\expInfo';

ee = [b.evts];
pp = [b.params];
ntr = b.completedTrials;
Fs  = 1/median(diff(t));

baseline.firstTrialExcluded = false;        % set default to false - will turn true for some old ChoiceWorld datasets

%% note parameters of the baseline
% what happened last, what happens next, whether quiescence was imposed
disp('getting baseline parameters...');
baseline.params = [];

if b.contRecording
    % figure out if feedback and stimulus offset happened at the same time
    posF = find(ee.feedbackValues(1:ntr)>0);
    negF = find(ee.feedbackValues(1:ntr)<0);
    posFendtimes = ee.feedbackTimes(posF) + pp.rewardDuration;
    negFendtimes = ee.feedbackTimes(negF) + pp.negFeedbackDuration;
    feedbackEndTimes = zeros(1,ntr);
    feedbackEndTimes(posF) = posFendtimes;
    feedbackEndTimes(negF) = negFendtimes;
    
    if abs(mean(ee.stimuliOffTimes - feedbackEndTimes)) < 0.01
        baseline.params.preBaseline.stimOffset  = true;
        baseline.params.preBaseline.feedbackEnd = true;
        xx = mean(ee.endTrialTimes - feedbackEndTimes);
        if abs(xx) < 0.01
            trialEndITI = false;
        else
            trialEndITI = true;
        end
    else
        if  mean(ee.stimuliOffTimes - feedbackEndTimes) > 0
            baseline.params.preBaseline.stimOffset  = true;
            baseline.params.preBaseline.feedbackEnd = false;
            xx = mean(ee.endTrialTimes - ee.stimuliOffTimes);
            if abs(xx) < 0.01
                trialEndITI = false;
            else
                trialEndITI = true;
            end
        else baseline.params.preBaseline.stimOffset = false;
            baseline.params.preBaseline.feedbackEnd = true;
            xx = mean(ee.endTrialTimes - feedbackEndTimes);
            if abs(xx) < 0.01
                trialEndITI = false;
            else
                trialEndITI = true;
            end
        end
    end
    baseline.params.preBaseline.recordingOnset      = false;
else
    baseline.params.preBaseline.stimOffset      = false;
    baseline.params.preBaseline.feedbackEnd     = false;
    baseline.params.preBaseline.recordingOnset  = true;
end

% check if there was an onset and go cue, check if one of them was zero
% for the purposes of baseline finding, only the cue with sound matters
if size(b.evts.cueOnTimes,1) > 1 
    if b.params.cueParams.ToneAmp(1) > 0
        switch b.params.cueType{2}
            case 'goCue'
                ee.cueOnTimes = b.evts.cueOnTimes(2,1:ntr);
            case 'preStimCue'
                ee.cueOnTimes = b.evts.cueOnTimes(1,1:ntr);
            case 'preStimCue+goCue'
                ee.cueOnTimes = b.evts.cueOnTimes(1,1:ntr);
        end
    end
end

% figure out what happened at end of baseline period, stimulus or cue or both?
if mean(([ee.stimuliOnTimes(1:ntr)]-[ee.cueOnTimes(1:ntr)])) < 0          % stimulus preceded cue
    baseline.params.baselineEnd.stimulusOnset   = true;
    baseline.params.baselineEnd.cueOnset        = false;
elseif mean(([ee.stimuliOnTimes(1:ntr)]-[ee.cueOnTimes(1:ntr)])) == 0     % stimulus and cue came on simultaneously
    baseline.params.baselineEnd.stimulusOnset   = true;
    baseline.params.baselineEnd.cueOnset        = true;
else
    baseline.params.baselineEnd.stimulusOnset   = false;
    baseline.params.baselineEnd.cueOnset        = true;
end


% figure out whether quiescence was imposed
if mean(pp.preStimQP) > 0
    baseline.params.baselineEnd.stimulusQuiesc  = true;
else baseline.params.baselineEnd.stimulusQuiesc = false;
end
if mean(pp.stimQP) > 0
    baseline.params.baselineEnd.quiescImposed   = true;
else baseline.params.baselineEnd.quiescImposed  = false;
end
    
%% find baseline onset & offset for each trial
disp('finding baseline onsets and offsets...')

baselineOnOffTimes = zeros(ntr,2);       % first column is baseline start, second column is baseline end
if b.contRecording
    baselineOnOffTimes(1,1)  = ee.newTrialTimes(1);     % first trial will always start here
    if trialEndITI      % use last feedback/stim offset as beginning
        if baseline.params.preBaseline.stimOffset
            baselineOnOffTimes(2:end,1) = ee.stimuliOffTimes(1:ntr-1) + 1/Fs;
        elseif baseline.params.preBaseline.feedbackEnd
            baselineOnOffTimes(2:end,1) = feedbackEndTimes + 1/Fs;
        else error('preBaseline parameters not correctly defined');
        end
    else                % use trial end timing as beginning of baseline
        baselineOnOffTimes(2:end,1) = ee.endTrialTimes(1:ntr-1);
    end
    
    if baseline.params.baselineEnd.stimulusOnset
        baselineOnOffTimes(:,2) = ee.stimuliOnTimes(1:ntr) - 1/Fs;
    elseif baseline.params.baselineEnd.cueOnset
        baselineOnOffTimes(:,2) = ee.cueOnTimes(1:ntr) - 1/Fs;
    else error('baselineEnd parameters not correctly defined');
    end
    
else
    % load timeline
    load(fullfile(BehDir,b.animal,b.iseries,b.iexp,strcat(b.iseries,'_',b.iexp,'_',b.animal,'_Timeline.mat')));
    tlOffset = Timeline.currSysTimeTimelineOffset;
    
    if baseline.params.preBaseline.recordingOnset
        % find baseline onset by finding cam2 onset after each break in imaging
        DAQts = Timeline.rawDAQTimestamps;
        cam2 = Timeline.rawDAQData(:, Timeline.hw.inputs(...
            strcmp({Timeline.hw.inputs.name},'cam2')).arrayColumn);
        
        % find camera frame timestamps
        cam2diff        = zeros(1,length(cam2));
        cam2diff(2:end) = diff(cam2);
        cam2OnsetsAll   = find(cam2diff>1); %finds the index number within cam2 in which the value jumps up, which corresponds to onset of imaging
        cam2FrameTimes  = DAQts(cam2OnsetsAll);
        clear cam2diff
        
        c2diff      = zeros(1,length(cam2FrameTimes));
        c2diff(1)   = (5);  % arbitrarily set this to higher than 1, as the first value in cam2FrameTimes is the onset of the first trial
        c2diff(2:end) = diff(cam2FrameTimes);
        cam2TrialOnsets = find(c2diff>1); %finds the indices of the trial onsets within cam2FrameTimes
        cam2TrialOnsetFrameTimes = cam2FrameTimes(cam2TrialOnsets);
        clear c2diff
        
        switch b.expType
            case 'ChoiceWorld'
                if b.iseries(4) == '5'        % quick hack... most likely this dataset is one of Elina's first ChoiceWorld datasets
                    % then the first trial wasn't recorded
                    baselineOnOffTimes(1) = b.evts.newTrialTimes(1);
                    baselineOnOffTimes(2:end,1) = cam2TrialOnsetFrameTimes(1:ntr-1)';
                end
            case 'signals'
                if length(cam2TrialOnsetFrameTimes) == ntr+1             % last trial incomplete but imaging started
                    baselineOnOffTimes(1:end,1) = cam2TrialOnsetFrameTimes(1:end-1)';
                else
                    baselineOnOffTimes(1:end,1) = cam2TrialOnsetFrameTimes';
                end
        end
    else
        error('preBaseline parameters not correctly defined');
    end
    
    if baseline.params.baselineEnd.stimulusOnset
        baselineOnOffTimes(:,2) = ee.stimuliOnTimes(1:ntr) - 1/Fs;
    elseif baseline.params.baselineEnd.cueOnset
        baselineOnOffTimes(:,2) = ee.cueOnTimes(1:ntr) - 1/Fs;
    else error('baselineEnd parameters not correctly defined');
    end
    
    % error('you havent written the case for non-continuous datasets yet!');
end

if baseline.params.preBaseline.recordingOnset    % if recording wasn't continuous, chop off first 200ms at imaging onset as onset causes artefacts
    baselineOnOffTimes(:,1) = baselineOnOffTimes(:,1) + 0.2;
    disp('Non-continuous imaging dataset, chopping off first 200ms of baseline to get rid of onset artefacts');
end

xx = diff(baselineOnOffTimes');
if min(xx)<0
    error('one of the baseline start times is after the baseline end time!');
end

%% find corresponding frames in t/V
% gets indices of frames in V, timepoints in t
bl_Vframes  = interp1(t, 1:numel(t), baselineOnOffTimes, 'nearest');
while isnan(bl_Vframes(1)) || isnan(bl_Vframes(1,2))              % if behavioural block started before recording
    if t(1) > baselineOnOffTimes(1,2)       % first baseline wasn't imaged, therefore exclude first trial
        ntr = ntr-1;
        baselineOnOffTimes = baselineOnOffTimes(2:end,:);
        bl_Vframes         = bl_Vframes(2:end,:);
        baseline.firstTrialExcluded = true;
    else
        bl_Vframes(1)       = 1;
        baselineOnOffTimes(1)    = t(1); % to start the baseline when the actual imaging started
    end
end
while isnan(bl_Vframes(end))         % if behavioural block kept going even though imaging ended
    bl_Vframes          = bl_Vframes(1:end-1,:);
    baselineOnOffTimes  = baselineOnOffTimes(1:end-1,:);
    ntr = ntr-1;
end

xx = diff(bl_Vframes');
maxBLfr = max(xx)+1; clear xx     % this gives the maximum number of frames within a baseline

baseline.OnOffTimes = baselineOnOffTimes;
baseline.V_OnOffFramenumbers = bl_Vframes;

%% 
% find corresponding indices of time points in wheel trace
bl_wheelInd = interp1(b.wheel.Times, 1:numel(b.wheel.Times), baselineOnOffTimes, 'nearest');
if isnan(bl_wheelInd(1))
    bl_wheelInd(1) = 1;
end

%% find corresponding frames in pupil/face movie
if doPupil
    [fPupilSize, eFT, doPupil] = loadPupil(b);
end

if doPupil

    pupilNaNs = find(isnan(fPupilSize));
    isPupilNaN = zeros(length(eFT),1);
    isPupilNaN(pupilNaNs) = 1;
    pupilNan_times = eFT(pupilNaNs);
    
    fPS_finite = fPupilSize(isfinite(fPupilSize));
    eFT_finite = eFT(isfinite(fPupilSize));
    fPS_interp = interp1(eFT_finite,fPS_finite,eFT);
    
    pupil_t = interp1(eFT, fPS_interp, t, 'linear');
    % get indices within t that we need to set to NaN in pupil_t
    isPupilNaN_int = interp1(eFT,isPupilNaN,t,'linear')';
    isPupilNaN_int(isnan(isPupilNaN_int)) = 0;      % if t starts earlier than eFT, there will be NaN values which we want to set to zero
    isPupilNaN_int = logical(isPupilNaN_int);
    pupil_t(isPupilNaN_int) = NaN;
    
%     % to double check that the fit is good, plot figure
%     figure;
%     plot(eFT, fPupilSize, 'b.'); hold on;
%     plot(t,pupil_t,'yo');
%     xlabel('time (s)');
%     ylabel('pupil area');
%     legend({'pupil','pupil interpolated to imaging time'});
%     title([b.animal ' ' b.iexp ' ' b.iseries]);
       
end
%%
% interpolate wheel and pupil to match timestamps in t
disp('finding baseline wheel and pupil values...');

baselineWheel = NaN(ntr,maxBLfr);   % each row is a trial
baselinePupil = NaN(ntr,maxBLfr);
baselineTimes = NaN(ntr,maxBLfr);

wheel_t = interp1(b.wheel.Times, b.wheel.Values, t);
for iit = 1:ntr
    thisT   = t(bl_Vframes(iit,1):bl_Vframes(iit,2));
    baselineTimes(iit,1:length(thisT)) = thisT;
    
    % find wheel trajectory during baseline time of trial iit
    thisWh  = wheel_t(bl_Vframes(iit,1):bl_Vframes(iit,2));
    baselineWheel(iit,(maxBLfr-length(thisT))+1:end) = thisWh - thisWh(end);        % set each baseline to start at 0 position for wheel
    
    if doPupil
        % find pupil dynamics
        % interpolate entire trace of pupil to match imaging timepoints
        % then just use indices of t timepoints to get pupil
        thisPp = pupil_t(bl_Vframes(iit,1):bl_Vframes(iit,2));
        baselinePupil(iit,(maxBLfr-length(thisT))+1:end) = thisPp;      % to have them all aligned to end of baseline, as beginning of baseline is more various
    end
    
end

baseline.Vtimepoints    = baselineTimes;
baseline.wheelValues    = baselineWheel;
baseline.pupilTraces    = baselinePupil;

disp('done');

end