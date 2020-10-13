% script that otherwise does the same as behaviourPerPower_imaging, but
% uses electrophysiology data instead

% written by Elina Jacobs, UCL Cortexlab

% set directories here
thisDir = myDirectoryWithTheData;

% requires Chronux toolbox

stimName = 'Visual_ephys';
frb = [3 6];    % frequency range of interest
doPupil = false;

minBL = 0.7;    % minimum baseline length

% load the list with the ephys experiments we want to look at
load(fullfile(myDirectoryWithTheData,'ephys_expList_forPowerAnalysis.mat'));
AL = expList_forPowerAnalysis; clear expList_forPowerAnalysis

nExps = size(AL,2);

nPerc   = 5;       % how many percentiles to divide power into
percentiles = linspace(0,100,nPerc+1);
% initialise matrices:
% for comparison of all trials that had a stimulus
anyC_percGO     = zeros(length(AL),nPerc)*NaN;
anyC_percNOGO   = zeros(length(AL),nPerc)*NaN;
% for comparison of all trials that had a stimulus, excluding the ones with equal contrasts
anyC_percGO_noEC   = zeros(length(AL),nPerc)*NaN;
anyC_percNOGO_noEC = zeros(length(AL),nPerc)*NaN;
anyC_percCorr_noEC = zeros(length(AL),nPerc)*NaN;
anyC_percInc_noEC  = zeros(length(AL),nPerc)*NaN;
% for high contrast trials (with contrasts on both sides as well)
highC_percGO    = zeros(length(AL),nPerc)*NaN;
highC_percNOGO  = zeros(length(AL),nPerc)*NaN;
highC_percCorr  = zeros(length(AL),nPerc)*NaN;
highC_percInc   = zeros(length(AL),nPerc)*NaN;
% for low contrast trials (with contrasts on both sides as well)
lowC_percGO    = zeros(length(AL),nPerc)*NaN;
lowC_percNOGO  = zeros(length(AL),nPerc)*NaN;
lowC_percCorr  = zeros(length(AL),nPerc)*NaN;
lowC_percInc   = zeros(length(AL),nPerc)*NaN;
% for zero contrast trials
zeroC_percCNG   = zeros(length(AL),nPerc)*NaN;
zeroC_percFA    = zeros(length(AL),nPerc)*NaN;

%% loop through experiments
for iie=1:nExps
    
    V1data = false;
    
    Exps = AL{iie};
    if ~isfield(Exps,'animal')
        Exps.animal = Exps.mousename;
    end
    if ~isfield(Exps,'iseries')
        Exps.iseries = Exps.series;
    end
    if ~isfield(Exps,'iexp')
        Exps.iexp = Exps.exp;
    end
    expRef = strcat(Exps.iseries(1:4),'-',Exps.iseries(5:6),'-',Exps.iseries(7:8),...
        '_',Exps.iexp,'_',Exps.animal);
    
    
    % figure out whether V1 was recorded from in this experiment
    ROIs = {};  % initialise the ROIs for each experiment anew
    thisEphysDir = fullfile(ephysDir,Exps.animal,...
        strcat(Exps.iseries(1:4),'-',Exps.iseries(5:6),'-',Exps.iseries(7:8)));
    ff = dir(thisEphysDir);
    folders = [ff(3:end).isdir];
    for xx = 1:length(folders)
        if folders(xx)
            folderName = ff(2+xx).name;
            if ~isempty(strfind(folderName,'ephys'))
                if isequal(strrep(folderName,'ephys_',''),'V1')
                    V1data = true;
                end
            end
        end
    end
    
    if V1data
        %% load data
        
        % load behavioural data
        [b] = generateGenBlock(expRef, Exps);
        ntr = b.completedTrials;
        
        % get behavioural variables
        [rGO, rNOGO, rCorr, rInc, rCNG, rGOc, rIncZC] = ...
            getChoiceNeglectInds(b,ntr);
        
        if isfield(b,'excludeFirstTrial')
            if b.excludeFirstTrial
                stim = b.stimuli(2:ntr+1,:);
            else
                stim  = b.stimuli(1:ntr,:);
            end
        else
            stim  = b.stimuli(1:ntr,:);
        end
        
        x1 = ismember(stim,[0 0],'rows');
        zeroCtrials = find(x1);
        clear x1
        
        equalCtrials    = find(diff(stim')==0);
        equalCNZtrials  = setdiff(equalCtrials,zeroCtrials);
        
        bothCtrials     = intersect(find(stim(:,1)>0),find(stim(:,2)>0));
        
        [highCtrials,~] = find(stim>=0.5); highCtrials = sort(highCtrials);
        highCtrialsNE   = setdiff(highCtrials,equalCNZtrials);
        highCtrialsNCC  = setdiff(highCtrials,bothCtrials);
        
        [lowCtrials_any,~]  = find(stim<0.5 & stim~=0);
        lowCtrials      = setdiff(lowCtrials_any,highCtrials);
        lowCtrialsNE    = setdiff(lowCtrials,equalCNZtrials);
        lowCtrialsNCC   = setdiff(lowCtrials,bothCtrials);
        
        anyCtrials      = setdiff(1:size(stim,1),zeroCtrials);
        anyCtrialsNOEC  = setdiff(1:size(stim,1),equalCtrials);
        
        
        if ~isempty(intersect(zeroCtrials,lowCtrials)) && ~isempty(intersect(zeroCtrials,highCtrials)) ...
                && ~isempty(intersect(zeroCtrials,anyCtrials))
            error('something is wrong with the contrast indexing');
        end
        
        %% load ephys data to compute power per Trial
        thisEphysDir = fullfile(myDirectoryWithTheData, b.animal,b.iseries,strcat('ephys_V1'),'sorting');
        [spikeTimes, spikeAmps, spikeDepths, spikeSites ] = ksDriftmap(thisEphysDir);
        
        % load manually selected channel numbers for V1
        load(fullfile(myDirectoryWithTheData,'manually_selected_V1_cortical_channel_borders_for_ephys_data',...
            strcat('V1_cortexBorders_',b.animal,'_',b.iseries,'_',b.iexp,'.mat')));
        
        % get indices of spikes that are in the cortical area we're interested in
        theseSpikeInds = ...
            find(spikeSites>=channelBoundaries(1)&spikeSites<=channelBoundaries(2));
        theseSpikeTimes = spikeTimes(theseSpikeInds);
        
        theseSpikeAmps  = spikeAmps(theseSpikeInds);
        % look at the spike data at same resultion as imaging, which is 35Hz
        imFs = 35;
        gaussWinLength = 6;
        
        MUA = histc(theseSpikeTimes, 0:(1/imFs):max(theseSpikeTimes));  % generate MUA
        gaussWin = gausswin(gaussWinLength); gaussWin = gaussWin/sum(gaussWin);      % sum of gaussian window needs to 1 so that it doesn't increase/decrease firing rate artificially
        smoothMUA= conv(MUA, gaussWin);
        smoothMUA = smoothMUA((gaussWinLength/2):end-(gaussWinLength/2));
        MUAt = linspace((1/imFs),max(theseSpikeTimes),length(MUA));     % generate time vector for MUA
        
        [baseline, ~] = get_baseline(b, MUAt, false);
        [baseline]  = find_noMvmtBaseline(baseline, ntr, imFs);
        shortbl     = find(diff(baseline.noMvmtFrames')<minBL*imFs);
        
        %% compute power per Trial
        powerPerTrial = zeros(ntr,1)*NaN;
        for iit = 1:ntr
            if isempty(shortbl)
                thisMUA = smoothMUA(baseline.noMvmtFrames(iit,1):baseline.noMvmtFrames(iit,2));
                msMUA   = thisMUA - mean(thisMUA);  % subtract mean from the MUA trace
                powerPerTrial(iit) = bandpower(msMUA,imFs,frb);
            else
                % for some reason, bandpower still computes the power even if
                % the baseline segment is too short - therefore 'manually
                % avoid' looking at those trials
                if ~ismember(shortbl,iit)
                    thisMUA = smoothMUA(baseline.noMvmtFrames(iit,1):baseline.noMvmtFrames(iit,2));
                    msMUA   = thisMUA - mean(thisMUA);  % subtract mean from the MUA trace
                    powerPerTrial(iit) = bandpower(msMUA,imFs,frb);
                end
            end
        end
        
      
        %% compute power per percentile
        powerPercentiles = prctile(powerPerTrial,percentiles);
        
        for prc = 1:nPerc
            thisPrcTrials = find(powerPerTrial>=powerPercentiles(prc) & ...
                powerPerTrial<powerPercentiles(prc+1));
            
            for nc = 1:7    % loop through the 4 stimulus conditions
                switch nc
                    case 1
                        theseTrials = intersect(thisPrcTrials,anyCtrials);
                        if ~isempty(theseTrials)
                            if length(intersect(theseTrials,rGO))+length(intersect(theseTrials,rNOGO)) ...
                                    ~= length(theseTrials)
                                error('something does not add up');
                            end
                            anyC_percGO(iie,prc) = ...
                                length(intersect(theseTrials,rGO))/length(theseTrials);
                            anyC_percNOGO(iie,prc) = ...
                                length(intersect(theseTrials,rNOGO))/length(theseTrials);
                        end
                        
                    case 2
                        theseTrials = intersect(thisPrcTrials,anyCtrialsNOEC);
                        if ~isempty(theseTrials)
                            if length(intersect(theseTrials,rCorr)) + length(intersect(theseTrials,rInc)) ...
                                    + length(intersect(theseTrials,rNOGO)) ~= length(theseTrials)
                                error('something does not add up');
                            end
                            anyC_percGO_noEC(iie,prc)   = ...
                                length(intersect(theseTrials,rGO))/length(theseTrials);
                            anyC_percNOGO_noEC(iie,prc)   = ...
                                length(intersect(theseTrials,rNOGO))/length(theseTrials);
                            
                            theseTrialsCorr = intersect(theseTrials,rCorr);
                            theseTrialsInc  = intersect(theseTrials,rInc);
                            theseTrialsSum  = length(theseTrialsCorr)+length(theseTrialsInc);
                            if theseTrialsSum > 0
                                anyC_percCorr_noEC(iie,prc) = ...
                                    length(theseTrialsCorr)/theseTrialsSum;
                                anyC_percInc_noEC(iie,prc) = ...
                                    length(theseTrialsInc)/theseTrialsSum;
                            end
                        end
                        
                    case 3
                        theseTrials = intersect(thisPrcTrials,highCtrialsNE);
                        if ~isempty(theseTrials)
                            if length(intersect(theseTrials,rCorr))+length(intersect(theseTrials,rInc)) ...
                                    + length(intersect(theseTrials,rNOGO)) ~= length(theseTrials)
                                error('something does not add up');
                            end
                            highC_percGO(iie,prc) = ...
                                length(intersect(theseTrials,rCorr))/length(theseTrials);
                            highC_percNOGO(iie,prc) = ...
                                length(intersect(theseTrials,rNOGO))/length(theseTrials);
                            
                            theseTrialsCorr = intersect(theseTrials,rCorr);
                            theseTrialsInc  = intersect(theseTrials,rInc);
                            theseTrialsSum  = length(theseTrialsCorr)+length(theseTrialsInc);
                            if theseTrialsSum > 0
                                highC_percCorr(iie,prc) = ...
                                    length(theseTrialsCorr)/theseTrialsSum;
                                highC_percInc(iie,prc) = ...
                                    length(theseTrialsInc)/theseTrialsSum;
                            end
                        end
                        
                    case 4
                        theseTrials = intersect(thisPrcTrials,lowCtrialsNE);
                        if ~isempty(theseTrials)
                            if length(intersect(theseTrials,rCorr))+length(intersect(theseTrials,rInc)) ...
                                    + length(intersect(theseTrials,rNOGO)) ~= length(theseTrials)
                                error('something does not add up');
                            end
                            lowC_percGO(iie,prc) = ...
                                length(intersect(theseTrials,rCorr))/length(theseTrials);
                            lowC_percNOGO(iie,prc) = ...
                                length(intersect(theseTrials,rNOGO))/length(theseTrials);
                            
                            theseTrialsCorr = intersect(theseTrials,rCorr);
                            theseTrialsInc  = intersect(theseTrials,rInc);
                            theseTrialsSum  = length(theseTrialsCorr)+length(theseTrialsInc);
                            if theseTrialsSum > 0
                                lowC_percCorr(iie,prc) = ...
                                    length(theseTrialsCorr)/theseTrialsSum;
                                lowC_percInc(iie,prc) = ...
                                    length(theseTrialsInc)/theseTrialsSum;
                            end
                        end
                        
                    case 5
                        theseTrials = intersect(thisPrcTrials,zeroCtrials);
                        if ~isempty(theseTrials)
                            if length(intersect(theseTrials,rCNG))+length(intersect(theseTrials,rIncZC)) ...
                                    ~= length(theseTrials)
                                error('something does not add up');
                            end
                            zeroC_percCNG(iie,prc) = ...
                                length(intersect(theseTrials,rCNG))/length(theseTrials);
                            zeroC_percFA(iie,prc) = ...
                                length(intersect(theseTrials,rIncZC))/length(theseTrials);
                        end
                end
                
            end
            
        end
    end
end


