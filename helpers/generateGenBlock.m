function [b] = generateGenBlock(varargin)
% script to generate 'general purpose block' from both signals and choiceworld inputs
% align timepoints to Timeline, if it exists

expRef = varargin{1};
mn  = expRef(14:end);
is  = expRef(1:10);
ie  = expRef(12);       % if exp > 9, this hack returns incorrect mn, is and ie - therefore the two cases
% load block
BehDir = '\\zserver.cortexlab.net\Data\expInfo';
if nargin < 2
    blockName = fullfile(BehDir,mn,is,ie,strcat(expRef,'_Block.mat'));
    tlName    = fullfile(BehDir,mn,is,ie,strcat(expRef,'_Timeline.mat'));
elseif nargin < 3
    Exps = varargin{2};
    blockName = fullfile(BehDir, Exps.animal, is, Exps.iexp, strcat(expRef,'_Block.mat'));
    tlName    = fullfile(BehDir, Exps.animal, is, Exps.iexp, strcat(expRef,'_Timeline.mat'));
elseif nargin == 3
    % case for early auditory datasets
    Exps = varargin{2};
    OldAudBehDir = fullfile(BehDir,'AVOld');
    blockName = fullfile(OldAudBehDir,strcat(Exps.iseries,'_',Exps.animal,'_',Exps.iexp,'.mat'));
    tlName    = fullfile(OldAudBehDir,strcat(expRef,'_Timeline.mat'));
else error('Too many input arguments');
end

% initialise general purpose block
b           = [];
if exist('Exps')
    mn = Exps.animal;
    ie = Exps.iexp;
end
b.animal    = mn;
b.iseries   = is;
b.iexp      = ie;

b.evts      = [];
b.params    = [];
b.stimuli   = [];
b.wheel     = [];

b.completedTrials= [];
b.recording      = false;       % set defaults to false
b.contRecording  = false;       % set defaults to false
b.indicator      = [];
b.rigName        = [];
b.expType        = [];

if exist(blockName)
    load(blockName);
    
    % following are the same in both block types
    if nargin == 3
        b.expStartedTime = block.StartDateTime;         
        b.expEndedTime   = [];
    else
        b.expStartedTime = [block.experimentStartedTime];
        b.expEndedTime   = [block.experimentEndedTime];
    end
    
    if isfield(block,'expType')
        expType = [block.expType];
    else
        if isfield(block,'expDef')
            expType = 'signals';
        elseif nargin == 3
            expType = 'Old_Auditory';
        else
            error('unknown block type');
        end
    end
    
    switch expType
        
        case 'ChoiceWorld'
            
            b.expType = 'ChoiceWorld';
            
            cc   = [block.trial.condition];
            nctr = [block.numCompletedTrials];
            if block.numCompletedTrials > 0
                % fill events sub-structure
                b.evts.newTrialTimes    = [block.trial.trialStartedTime];
                cot = [block.trial.onsetToneSoundPlayedTime];
                if length(cot) > 1.5*nctr
                    if length(cot)/2 > nctr             % last incomplete trial included the onset cue
                        if length(cot) == nctr*2 + 1
                            cot(end+1) = NaN;
                            cot = reshape(cot,2,nctr+1);
                        else
                            cot = reshape(cot,2,nctr+1);
                        end
                    else cot = reshape(cot,2,nctr);
                    end
                    b.evts.cueOnTimes   = cot;      % first row preStimCue, second row goCue
                else b.evts.cueOnTimes  = [block.trial.onsetToneSoundPlayedTime];
                end
                b.evts.stimuliOnTimes   = [block.trial.stimulusCueStartedTime];
                b.evts.stimuliOffTimes  = [block.trial.stimulusCueEndedTime];
                b.evts.interactiveOnT   = [block.trial.interactiveStartedTime];
                b.evts.responseTimes    = [block.trial.responseMadeTime];
                % in choiceworld, response IDs are different, 3 = nogo response, 1 vis stim on left , 2 vis stim on right
                % turn them into -1 for left, 1 for right, 0 for nogo
                rV = [block.trial.responseMadeID];
                rV(rV==3) = 0; rV(rV==1) = -1; rV(rV==2) = 1;
                b.evts.responseValues   = rV;
                b.evts.feedbackTimes    = [block.trial.feedbackStartedTime];
                b.evts.feedbackValues   = [block.trial.feedbackType];
                b.evts.endTrialTimes    = [block.trial.trialEndedTime];
                b.evts.repeatNum        = [cc.repeatNum];
                
                % if timeline exists, correct the time points to that in Steinmetz
                % dataset
                success = false;
                switch block.rigName
                    case 'zgood'
                        if exist(tlName);
                            load(tlName);
                            tt = Timeline.rawDAQTimestamps;
                            pd = Timeline.rawDAQData(:, strcmp({Timeline.hw.inputs.name}, 'photoDiode'));
                            pdT = schmittTimes(tt, pd, [3 4]);
                            sw = block.stimWindowUpdateTimes;
                            if length(sw)<=length(pdT) && length(sw)>1
                                [~,corr,success,actualTimes] = findCorrection(pdT, sw, false);
                            end
                            if success
                                events = {'newTrialTimes', 'cueOnTimes', 'stimuliOnTimes', 'stimuliOffTimes', ...
                                    'interactiveOnT', 'responseTimes', 'feedbackTimes', 'endTrialTimes'};
                                for e = 1:length(events)
                                    b.evts.([events{e}]) = b.evts.([events{e}]) + corr(2);
                                end
                            end
                        end
                end
                
                % fill stimuli sub-structure
                ss = [cc.visCueContrast];
                b.stimuli               = ss';      % this was, left column corresponds to left stim, right col to right stim
                
                % fill parameters sub-structure
                b.params.ITI            = [block.parameters.interTrialDelay];
                b.params.preStimQP      = [block.parameters.preStimQuiescentPeriod];
                b.params.stimQP         = 0;    % this option does not usually exist in ChoiceWorld
                % CW doesn't have visual cue option, so just check the ToneAmp is not zero
                if isfield(block.parameters, 'onsetToneMaxAmp')     % then it's an old dataset
                    if block.parameters.onsetToneMaxAmp > 0
                        b.params.cueType{1}     = 'Tone';
                        if median([b.evts.cueOnTimes] - [b.evts.stimuliOnTimes]) < 0
                            b.params.cueType{2} = 'preStimCue';
                        elseif median([b.evts.cueOnTimes] - [b.evts.stimuliOnTimes]) > 0
                            b.params.cueType{2} = 'goCue';
                        else
                            b.params.cueType{2} = 'simultaneousCue&Stim';
                        end
                    else
                        b.params.cueType{1} = [];
                        b.params.cueType{2} = [];
                    end
                    b.params.cueParams.ToneAmp = [block.parameters.onsetToneMaxAmp];
                else
                    if block.parameters.onsetToneRelAmp > 0 || block.parameters.interactiveOnsetToneRelAmp > 0
                        b.params.cueType{1}     = 'Tone';
                        if block.parameters.onsetToneRelAmp > 0 && block.parameters.interactiveOnsetToneRelAmp == 0
                            b.params.cueType{2} = 'preStimCue';
                        elseif block.parameters.onsetToneRelAmp == 0 && block.parameters.interactiveOnsetToneRelAmp > 0
                            b.params.cueType{2} = 'goCue';
                        else b.params.cueType{2}= 'preStimCue+goCue';
                        end
                    else
                        b.params.cueType{1} = [];
                        b.params.cueType{2} = [];
                    end
                    b.params.cueParams.ToneAmp  = [block.parameters.onsetToneRelAmp,block.parameters.interactiveOnsetToneRelAmp]; 
                end
                b.params.cueDuration    = max([block.parameters.onsetToneDuration]);
                b.params.cueParams.Contrast = 0;
                b.params.stimType{1}    = 'Visual';
                ds = diff(ss); ng = find(ds==0);
                if isempty(ng)
                    b.params.stimType{2}      = '2AFC';
                else
                    ss_ng = ss(1,:); ss_ng = ss_ng(ng);
                    if max(ss_ng) > 0
                        b.params.stimType{2}  = 'ContrastDiscrimination';
                    else
                        fng = [b.evts.feedbackValues(ng)];
                        rng = rV(ng);           % find responses during 0 contrast trials
                        ngr = find(rng==0);     % find the trials with nogo response during 0 contrast trials
                        if isempty(ngr) || mean(fng(ngr))==-1 % ie giving nogo response is always wrong during a 0 contrast trial, which is the case for Daisuke
                            b.params.stimType{2} = '2AFC';
                        else
                            ss1 = ss(1,:);
                            ss2 = ss(2,:);
                            c1  = find(ss1>0); % finds all trials in which there was a contrast on right
                            c2  = find(ss2>0);
                            cc2s = intersect(c1,c2);
                            if isempty(cc2s)    % then there were never two contrasts simultaneously
                                b.params.stimType{2} = '2ANFC';
                            else
                                b.params.stimType{2}  = 'ContrastDiscrimination';
                            end
                        end
                    end
                end
                b.params.responseWindow     = [block.parameters.responseWindow];
                b.params.rewardDuration     = [block.parameters.positiveFeedbackPeriod];
                b.params.rewardSize         = [block.parameters.rewardVolume];
                b.params.negFeedbackDuration= [block.parameters.negFeedbackSoundDuration];
                b.params.negFeedbackAmp     = [block.parameters.negFeedbackSoundAmp];
                
                b.wheel.Values      = [block.inputSensorPositions] - block.inputSensorPositions(1);
                b.wheel.Times       = [block.inputSensorPositionTimes];
                
                if success
                    b.wheel.Times = b.wheel.Times + corr(2);
                    if b.wheel.Times(1) < 0
                        b.wheel.Times  = b.wheel.Times(2:end);
                        b.wheel.Values = b.wheel.Values(2:end);
                    end
                end
                
                switch block.rigName
                    case 'zooropa'
                        isx = is(is~='-');
                        if exist(fullfile('\\zserver.cortexlab.net\Data\GCAMP',mn,isx,ie))
                            b.recording = true;
                            b.indicator = 'GCAMP';
                        elseif exist(fullfile('\\zserver.cortexlab.net\Data\Subjects',mn,is,ie))
                            b.recording     = true;
                            b.contRecording = true;
                            b.indicator     = 'GCAMP';  % set as default GCAMP and use case below to set to VSFP if appropriate
                            if exist(fullfile('\\zserver.cortexlab.net\Data\Subjects',mn,is,'ops.mat'))
                                load(fullfile('\\zserver.cortexlab.net\Data\Subjects',mn,is,'ops.mat'));
                                switch ops.userName
                                    case 'daisuke'
                                        b.indicator = 'VSFP';
                                end
                            end
                        end
                        b.rigName = 'BigRig';
                    case 'zgood'    % Nick's imaging rig
                        if exist(fullfile('\\zserver.cortexlab.net\Data\Subjects',mn,is,ie))
                            b.recording     = true;
                            b.contRecording = true;
                            b.indicator     = 'GCAMP';
                            b.rigName       = 'zgood';
                        elseif exist(fullfile('\\zserver.cortexlab.net\Data\Subjects',mn,is))
                            % ephys case - EJ addition 2017/12/06
                            b.recording     = true;
                            b.contRecording = true;
                            b.indicator     = 'ephys';
                            b.rigName       = 'zgood';
                        end
                end
                
                % if ephys recording, align all timestamps to ephys!
                % added by EJ on 2017/12/07
                if ~isempty(b.indicator)
                    switch b.indicator
                        case 'ephys'
                            SetDefaultDirs;
                            dataDir = fullfile(server2Name,'Data','Subjects');
                            expDir = fullfile(dataDir,b.animal,b.iseries,'alf');
                            if exist(expDir)
                                % load the aligned timestamps from npy files
                                ep.stimuliOnTimes   = readNPY(fullfile(expDir,'cwStimOn.times.npy'));
                                ep.responseTimes    = readNPY(fullfile(expDir,'cwResponse.times.npy'));
                                ep.feedbackTimes    = readNPY(fullfile(expDir,'cwFeedback.times.npy'));
                                % figure out what the correction factor is, check that it is consistent,
                                % and apply it to all timepoints of interest
                                delays.stimOn   = ep.stimuliOnTimes(1:nctr) - b.evts.stimuliOnTimes(1:nctr)';
                                delays.response = ep.responseTimes(1:nctr) - b.evts.responseTimes(1:nctr)';
                                delays.feedback = ep.feedbackTimes(1:nctr) - b.evts.feedbackTimes(1:nctr)';
                                
                                delaySTDs    = structfun(@std,delays);
                                delayMedians = structfun(@median,delays);
                                if max(delaySTDs) < 0.05 && std(delayMedians) < 0.05
                                    epCF = mean(delayMedians);
                                elseif max(delaySTDs) < 0.1 && std(delayMedians) < 0.1
                                    epCF = mean(delayMedians);
                                    warning('only around 100ms precision in time alignments between behaviour and ephys!');
                                else
                                    error('check dataset - time alignments dont match');
                                end
                                
                                % apply correction
                                b.evts.newTrialTimes    = b.evts.newTrialTimes + epCF;
                                b.evts.cueOnTimes       = b.evts.cueOnTimes + epCF;
                                b.evts.stimuliOnTimes   = b.evts.stimuliOnTimes + epCF;
                                b.evts.stimuliOffTimes  = b.evts.stimuliOffTimes + epCF;
                                b.evts.interactiveOnT   = b.evts.interactiveOnT + epCF;
                                b.evts.responseTimes    = b.evts.responseTimes + epCF;
                                b.evts.feedbackTimes    = b.evts.feedbackTimes + epCF;
                                b.evts.endTrialTimes    = b.evts.endTrialTimes + epCF;
                                b.wheel.Times           = b.wheel.Times +epCF;
                                
                            else
                                error('ephys dataset that cannot be aligned yet because alf folder does not exist');
                            end
                    end
                end
            end
            
            
        case 'signals'
            b.expType = 'signals';
            
            if exist(tlName)
                load(tlName);
                tlOffset = Timeline.currSysTimeTimelineOffset;
            else
                tlOffset = 0;
                warning('Timeline not found - if this is a recorded dataset timestamps will be wrong');
            end
            
            nctr = length(block.events.endTrialTimes);
            % fill events sub-structure
            b.evts.newTrialTimes    = [block.events.newTrialTimes] - tlOffset;
            b.evts.stimuliOnTimes   = [block.events.stimuliOnTimes] - tlOffset;
            if isfield(block.events,'cueOnTimes')
                b.evts.cueOnTimes   = [block.events.cueOnTimes] - tlOffset;
            else
                if isfield(block.paramsValues,'onsetToneDelay') && isfield(block.paramsValues,'onsetStimDelay')
                    cd = mean([block.paramsValues.onsetStimDelay]-[block.paramsValues.onsetToneDelay]);
                    b.evts.cueOnTimes = b.evts.stimuliOnTimes - cd;
                elseif  std([block.paramsValues.interactiveDelay]) == 0 && mean([block.paramsValues.interactiveDelay]) == 0
                    % this is a hack - if the interactiveDelay is zero, than the onset of that corresponds to onset of cue
                    if mean([block.events.interactiveOnTimes] - [block.events.stimuliOnTimes]) > 0  % then this was a block where cue followed stimulus
                        b.evts.cueOnTimes = [block.events.interactiveOnTimes] - tlOffset;
                    else
                        b.evts.cueOnTimes = [];
                    end
                elseif isfield(block.paramsValues,'onsetCueDelay') && isfield(block.paramsValues,'onsetStimDelay')
                    if min([block.paramsValues.onsetStimDelay]) == max([block.paramsValues.onsetStimDelay])
                        % then the delay between cue and stimulus was constant
                        cd = mean([block.paramsValues.onsetStimDelay]-[block.paramsValues.onsetCueDelay]);
                        b.evts.cueOnTimes = b.evts.stimuliOnTimes - cd;
                    else
                        b.evts.cueOnTimes = [];
                    end
                else
                    b.evts.cueOnTimes = [];    % cueOnTimes can be determined from delays, but first timeline offset needs to be substracted from all event times
                end
            end
            if isfield(block.events,'stimuliOffTimes')
                b.evts.stimuliOffTimes  = [block.events.stimuliOffTimes]  - tlOffset;
            elseif isfield(block.events,'pipOffTimes')
                b.evts.stimuliOffTimes  = [block.events.pipOffTimes]  - tlOffset;
            else
                error('unknown block type');
            end
            b.evts.interactiveOnT   = [block.events.interactiveOnTimes]  - tlOffset;
            b.evts.responseTimes    = [block.events.responseTimes] - tlOffset;
            b.evts.responseValues   = [block.events.responseValues];
            b.evts.feedbackTimes    = [block.events.feedbackTimes] - tlOffset;
            b.evts.feedbackValues   = [block.events.feedbackValues];
            b.evts.endTrialTimes    = [block.events.endTrialTimes] - tlOffset;
            b.evts.repeatNum        = [block.events.repeatNumValues];
            
            % fill stimuli sub-structure
            if isfield(block.paramsValues,'pipAmplitude')
                toneAmp             = [block.paramsValues.pipAmplitude];
            else toneAmp            = [];
            end
            if isfield(block.paramsValues,'targetContrast')
                visContrast         = [block.paramsValues.targetContrast];
            elseif isfield(block.paramsValues,'contrast')
                visContrast         = [block.paramsValues.contrast];
            elseif ~isempty(strfind(block.expDef,'tonePipWorld'));  % if early auditory signals script
                visContrast         = [];
            else error('unknown or inexistent contrast parameter');
            end
            azimuth                 = [block.paramsValues.targetAzimuth];
            ssL = find(azimuth<0); ssR = find(azimuth>0);
            ss = zeros(length(azimuth),2);
            if isempty(toneAmp)                                 % only visual stimuli
                ss(ssL,1) = visContrast(ssL);
                ss(ssR,2) = visContrast(ssR);
                b.stimuli           = ss;
                b.params.stimType{1} = 'Visual';
            elseif isempty(visContrast)                         % only auditory stimuli
                ss(ssL,1) = toneAmp(ssL);
                ss(ssR,2) = toneAmp(ssR);
                b.stimuli = ss;
                b.params.stimType{1} = 'Audio';
            else
                if max(visContrast) > 0 && max(toneAmp) == 0        % only visual stimuli
                    ss(ssL,1) = visContrast(ssL);
                    ss(ssR,2) = visContrast(ssR);
                    b.stimuli           = ss;
                    b.params.stimType{1} = 'Visual';
                elseif max(visContrast) == 0 && max(toneAmp) > 0    % only auditory stimuli
                    ss(ssL,1) = toneAmp(ssL);
                    ss(ssR,2) = toneAmp(ssR);
                    b.stimuli = ss;
                    b.params.stimType{1} = 'Audio';
                else ssV = ss; ssA = ss;
                    ssV(ssL,1) = visContrast(ssL);
                    ssV(ssR,2) = visContrast(ssR);
                    ssA(ssL,1) = toneAmp(ssL);
                    ssA(ssR,2) = toneAmp(ssR);
                    b.stimuli.visContrasts = ssV;
                    b.stimuli.toneAmps     = ssA;
                    b.params.stimType{1} = 'AudioVisual';
                end
            end
            %             if max(visContrast) > 0 && max(toneAmp) == 0        % only visual stimuli
            %                 ss(ssL,1) = visContrast(ssL);
            %                 ss(ssR,2) = visContrast(ssR);
            %                 b.stimuli           = ss;
            %                 b.params.stimType{1} = 'Visual';
            %             elseif max(visContrast) == 0 && max(toneAmp) > 0    % only auditory stimuli
            %                 ss(ssL,1) = toneAmp(ssL);
            %                 ss(ssR,2) = toneAmp(ssR);
            %                 b.stimuli = ss;
            %                 b.params.stimType{1} = 'Audio';
            %             else ssV = ss; ssA = ss;
            %                 ssV(ssL,1) = visContrast(ssL);
            %                 ssV(ssR,2) = visContrast(ssR);
            %                 ssA(ssL,1) = toneAmp(ssL);
            %                 ssA(ssR,2) = toneAmp(ssR);
            %                 b.stimuli.visContrasts = ssV;
            %                 b.stimuli.toneAmps     = ssA;
            %                 b.params.stimType{1} = 'AudioVisual';
            %             end
            
            % fill parameters sub-structure
            if isfield(block.paramsValues,'interTrialDelay')
                b.params.ITI        = [block.paramsValues(1).interTrialDelay];
            else b.params.ITI       = 0;
            end
            if isfield(block.paramsValues,'preStimQuiescentPeriod')
                b.params.preStimQP  = [block.paramsValues(1).preStimQuiescentPeriod];
            else b.params.preStimQP = 0;
            end
            if isfield(block.paramsValues,'stimQuiescentPeriod')
                b.params.stimQP     = [block.paramsValues(1).stimQuiescentPeriod];
            else b.params.stimQP    = 0;
            end
            if max([block.paramsValues.onsetToneAmplitude]) > 0
                if isfield(block.paramsValues,'cueContrast') && max([block.paramsValues.cueContrast]) > 0
                    b.params.cueType{1}  = 'Tone+Contrast';         % both tone and contrast were used as cue
                else b.params.cueType{1} = 'Tone';                  % only tone was used as cue
                end
            else
                if isfield(block.paramsValues,'cueContrast') && max([block.paramsValues.cueContrast]) > 0
                    b.params.cueType{1}  = 'Contrast';              % only contrast was used as cue
                else b.params.cueType{1} = [];                      % no cue was used
                end
            end
            if isfield(block.paramsValues,'onsetStimDelay') && isfield(block.paramsValues,'onsetToneDelay')
                if median([block.paramsValues.onsetStimDelay] - [block.paramsValues.onsetToneDelay]) > 0
                    b.params.cueType{2}     = 'preStimCue';
                elseif median([block.paramsValues.onsetStimDelay] - [block.paramsValues.onsetToneDelay]) < 0
                    b.params.cueType{2}     = 'goCue';
                else b.params.cueType{2}    = 'simultaneousCue&Stim';
                end
            elseif isfield(block.paramsValues,'onsetStimDelay') && isfield(block.paramsValues,'onsetCueDelay')
                if median([block.paramsValues.onsetStimDelay] - [block.paramsValues.onsetCueDelay]) > 0
                    b.params.cueType{2}     = 'preStimCue';
                elseif median([block.paramsValues.onsetStimDelay] - [block.paramsValues.onsetCueDelay]) < 0
                    b.params.cueType{2}     = 'goCue';
                else b.params.cueType{2}    = 'simultaneousCue&Stim';
                end
            else
                if mean([b.evts.stimuliOnTimes(1:nctr)] - [b.evts.cueOnTimes(1:nctr)]) > 0
                    b.params.cueType{2}     = 'preStimCue';
                elseif mean([b.evts.stimuliOnTimes(1:nctr)] - [b.evts.cueOnTimes(1:nctr)]) < 0
                    b.params.cueType{2}     = 'goCue';
                else b.params.cueType{2}    = 'simultaneousCue&Stim';
                end
            end
            %         if isempty(b.evts.cueOnTimes)
            %             b.params.cueType{2}     = 'preStimCue';             %
            %         else
            b.params.cueParams.ToneAmp      = max([block.paramsValues.onsetToneAmplitude]);
            if isfield(block.paramsValues,'cueContrast')
                b.params.cueParams.Contrast = max([block.paramsValues.cueContrast]);
            else b.params.cueParams.Contrast = 0;
            end
            if max([block.paramsValues.onsetToneAmplitude]) > 0
                b.params.cueDuration        = max([block.paramsValues.onsetToneDuration]);
            else
                if ~isempty(b.evts.cueOnTimes)
                    b.params.cueDuration  = [b.evts.stimuliOffTimes] - [b.evts.cueOnTimes(1:nctr)];
                else
                    b.params.cueDuration       = NaN;      % make this NaN given that the visual stimulus is always a go cue and stays until response time
                end
            end
            if max(max(ss)) > 0
                ds = diff(ss'); ng = find(ds==0);
                if isempty(ng)
                    b.params.stimType{2}    = '2AFC';
                else b.params.stimType{2}   = '2ANFC';
                end
            else
                dsV = diff(ssV'); ngV = find(dsV==0);
                dsA = diff(ssA'); ngA = find(dsA==0);
                % I never used two simultaneous stimuli, so true nogo trials would have both zero contrasts and zero amplitude,
                % otherwise they are single modality AFC trials within an audiovisual block
                ds = ngV - ngA;  ng  = find(ds==0);
                if isempty(ng)
                    b.params.stimType{2}    = '2AFC';
                else b.params.stimType{2}   = '2ANFC';
                end
            end
            b.params.responseWindow         = median([block.paramsValues.responseWindow]);
            b.params.rewardDuration         = median([block.paramsValues.rewardDur]);
            b.params.rewardSize             = median([block.paramsValues.rewardSize]);
            b.params.negFeedbackDuration    = median([block.paramsValues.noiseBurstDur]);
            b.params.negFeedbackAmp         = median([block.paramsValues.noiseBurstAmp]);
            
            b.wheel.Values      = [block.inputs.wheelValues] - block.inputs.wheelValues(1);
            b.wheel.Times       = [block.inputs.wheelTimes]  - tlOffset;
            
            %             if block.rigName == 'zooropa'
            switch block.rigName
                case 'zooropa'
                    if isfield(block.paramsValues,'onsetAcquisitionDelay')
                        % if this field exists, it was a recording, but not continuous
                        b.recording     = true;
                        b.indicator     = 'GCAMP';
                        b.rigName       = 'BigRig';
                    else
                        if exist(fullfile('\\zserver.cortexlab.net\Data\Subjects',block.expRef(14:18),block.expRef(1:10),block.expRef(12)))
                            b.recording     = true;
                            b.contRecording = true;
                            b.indicator     = 'GCAMP';
                            b.rigName       = 'BigRig';
                        elseif exist(fullfile('\\zserver.cortexlab.net\Data\Subjects',block.expRef(14:18),strcat(block.expRef(1:10),'_',block.expRef(12)),block.expRef(12)))
                            b.recording     = true;
                            b.contRecording = true;
                            b.indicator     = 'GCAMP';
                            b.rigName       = 'BigRig';
                        end
                    end
            end
            
            
        case 'Old_Auditory'
            b.expType = 'Old_Auditory';
            
            nctr = length([block.Trials.FeedbackStartTime]);
            % fill events sub-structure
            b.evts.newTrialTimes    = [block.Trials.IntermissionStartTime];
            b.evts.cueOnTimes       = [block.Trials.StimStartTime];      % there was no cue, so use stimOnset as cue Onset as otherwise code breaks, and specify in parameters that there was no actual cue
            b.evts.stimuliOnTimes   = [block.Trials.StimStartTime];
            b.evts.stimuliOffTimes  = [block.Trials.StimEndTime];
            b.evts.interactiveOnT   = [block.Trials.InteractiveStartTime];
            b.evts.responseTimes    = [block.Trials.InteractiveEndTime];
            b.evts.responseValues   = [block.Trials.ResponseDir];   % -1 for left/high stim, 1 for right/low stim
            b.evts.feedbackTimes    = [block.Trials.FeedbackStartTime];
            responseRewarded = [block.Trials.ResponseRewarded];
            responseRewarded(responseRewarded==0) = -1;
            b.evts.feedbackValues   = responseRewarded;
            b.evts.endTrialTimes    = [block.Trials.FeedbackEndTime];
            repeat = double([block.Trials.RepeatOfLast]);
            repeatNum = zeros(1,nctr);
            for iit = 1:nctr
                if repeat(iit) == 0
                    repeatNum(iit) = 1;
                else
                    if repeat(iit) == 1 && repeat(iit-1) == 0
                        repeatNum(iit) = 2;
                    else
                        repeatNum(iit) = repeatNum(iit-1)+1;
                    end
                end
            end
            b.evts.repeatNum        = repeatNum;
            
            % fill in stimulus sub structure
            rewardTargets   = reshape([block.Trials.RewardTargets],2,nctr)';
            ss              = bsxfun(@times,rewardTargets,[block.Trials.ToneAmplitude]');
            b.stimuli       = ss;
            
            % fill parameters sub-structure
            b.params.ITI       = [block.IntermissionMinDuration block.IntermissionMaxDuration];
            b.params.preStimQP = 0;     % I think? there is a separate 'preStimStart' and 'preStimEnd' time but not sure if that's quiescent...
            b.params.stimQP    = 0;
            b.params.cueType{1} = 'None'; % there was no cue
            b.params.cueParams.ToneAmp      = 0;
            b.params.cueParams.Contrast     = 0;
            b.params.cueDuration            = NaN;      % make this NaN given that the visual stimulus is always a go cue and stays until response time
            b.params.stimType{1}    = 'Audio';
            b.params.stimType{2}    = '2AFC';       % always 2AFC in old version
            
            if max([block.Trials.ResponseTimeout])
                error('check the response window length');
            else
                b.params.responseWindow = Inf;
            end
            b.params.rewardDuration         = block.PosFeedbackDuration;;
            b.params.rewardSize             = block.MainRewardSize;
            b.params.negFeedbackDuration    = block.NegFeedbackDuration;
            b.params.negFeedbackAmp         = block.AudNoiseAmplitude;
            
            % I believe this is the wheel?
            b.wheel.Values      = block.HelmAngles;
            b.wheel.Times       = block.HelmAngleTimes;
            
    end
    
    b.completedTrials       = nctr;
    b.excludeFirstTrial     = false;        % set it as false for a default, this can be changed after the baseline function?
    
else
    warning('block does not exist - returning empty structure');
end
end