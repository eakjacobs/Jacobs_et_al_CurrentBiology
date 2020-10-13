% script similar to summary_powerSpectrumRatio_imaging but for electrophysiology experiments
% generates powerspectrum ratio in visual cortex for every
% dataset; then plots average with SEM of that
%
% written by Elina Jacobs, UCL Cortexlab

% this script uses code from the spikes Github repository: 
% https://github.com/cortex-lab/spikes

% set directories here
thisDir = myDirectoryWithTheData;
% add path to Chronux toolbox

stimName = 'Visual_ephys';
frb = [3 6];    % frequency range of interest
doPupil = false;

minBL = 0.7;    % minimum baseline length
ec = 0;         % initialise index counter for summary structure

%%
% load the list with the ephys experiments we want to look at
load(fullfile(myDirectoryWithTheData,'ephys_expList_forPowerAnalysis_looserCriteria.mat'));
al = expList_forPowerAnalysis; clear expList_forPowerAnalysis

nExps = size(al,2);

% loop through experiments
for ee=1:nExps
    
    V1data = false;
    
    Exps = al{ee};
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
    
    
    % check whether V1 was recorded from in this experiment
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
        ec = ec+1;          % index counter for allRatios
        
        % load behavioural data
        [b] = generateGenBlock(expRef, Exps);
        ntr = b.completedTrials;
        [~, rNOGO, ~, ~, ~, rGOc, ~] = getChoiceNeglectInds(b,ntr,stimName);
        
        % load spike data
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
        
        % get baseline times
        [baseline, doPupil] = get_baseline(b, MUAt, doPupil);
        if baseline.firstTrialExcluded
            b.excludeFirstTrial = true;
        end
        [baseline]  = find_noMvmtBaseline(baseline, ntr, imFs);
        shortbl     = find(diff(baseline.noMvmtFrames')<minBL*imFs);
        
        %% compute powerspectra per trial
        
        params.tapers = [3 5];
        params.Fs     = round(1/median(diff(MUAt)));
        winSize       = minBL; % in seconds
        
        for iit = 1:ntr
            if isempty(shortbl)
                thisMUA = smoothMUA(baseline.noMvmtFrames(iit,1):baseline.noMvmtFrames(iit,2));
                msMUA   = thisMUA - mean(thisMUA);  % subtract mean from the MUA trace
                [p, f] = mtspectrumsegc(msMUA, winSize, params);
                allPowerSpectra(iit,:) = p';
            else
                % for some reason, bandpower still computes the power even if
                % the baseline segment is too short - therefore 'manually
                % avoid' looking at those trials
                if ~ismember(shortbl,iit)
                    thisMUA = smoothMUA(baseline.noMvmtFrames(iit,1):baseline.noMvmtFrames(iit,2));
                    msMUA   = thisMUA - mean(thisMUA);  % subtract mean from the MUA trace
                    [p, f] = mtspectrumsegc(msMUA, winSize, params);
                    allPowerSpectra(iit,:) = p';
                end
            end
        end
        
        allPowerSpectra(allPowerSpectra==0) = NaN;
        if sum(ismember(shortbl,ntr))
            allPowerSpectra(ntr,:) = NaN;
        end
        
        save(fullfile(myDirectoryWithTheData,...
            strcat(b.animal,'_',b.iseries,'_',b.iexp,'_Powerspectra_noMvmt_perTrial.mat')),'allPowerSpectra','f');
        
        %% fill in summary structure
        [~,maxInd] = min(abs(f-10));
        
        PSgo    = allPowerSpectra(rGOc,:);
        PSnogo  = allPowerSpectra(rNOGO,:);
        
        meanPSgo    = nanmean(PSgo,1);
        meanPSnogo  = nanmean(PSnogo,1);
        
        PSratio     = meanPSgo(1:maxInd) ./ meanPSnogo(1:maxInd);
        fss = linspace(f(1),f(maxInd),19);
        PSratio_ss = interp1(f(1:maxInd),PSratio,fss);
        
        allRatios(ec,:) = PSratio_ss;
        
        clear allPowerSpectra allPS PSgo PSnogo meanPSgo meanPSnogo PSratio PSratio_ss b t
        
        save(fullfile(myDirectoryWithTheData,...
            strcat(stimName,'_PowerRatio_noMvmt_notdB_VisualCortex.mat')),'allRatios','fss');
        
    end
    
end

%%
figure;
% plot([0 10],[1 1],'k:'); hold on;
shadedErrorBar(fss,nanmean(allRatios,1),nansem(allRatios,[],1),{'b-','markerfacecolor',[0.2 0.2 0.2]});
% hold on;
% plot(fss,allRatios(44,:),'b');
xlim([1 9]); 
xlabel('Frequency');
xticks([3 6 9]);
ylim([0.78 1.18]); 
yticks([0.8 0.9 1]);
ylabel('Choice / Neglect Power');
title('Visual Cortex Power Ratio');
box off;
set(gca, 'FontSize', 24);
axis square;

x0=200;
y0=200;
width = 400;
height= 400;
set(gca,'units','points','position',[x0,y0,width,height]);
