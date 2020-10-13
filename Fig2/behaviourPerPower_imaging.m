% script that bins power & calculates
% - percent correct choice
% - percent incorrect choice
% - percent neglect
% - percent incorrect go
% - percent correct nogo

% written by Elina Jacobs, UCL Cortexlab

clear all

% set directories here
thisDir = myDirectoryWithTheData;

stimName = 'VisualALL';

load(fullfile(myDirectoryWithTheData,...
    strcat(stimName,'_','expList_forPowerAnalysis.mat')));

AL = expList_forPowerAnalysis; clear expList_forPowerAnalysis

%%
frb     = [3,6];
nPerc   = 5;       % how many percentiles to divide power into
percentiles = linspace(0,100,nPerc+1);
switch stimName
    case 'Visual2ANFC'
        nROIs = 5;
    case 'VisualALL'
        nROIs = 5;
    case 'Auditory'
        nROIs = 4;
    case 'AudioVisual'
        nROIs = 4;
end

roiNames{1} = 'Visual';
roiNames{2} = 'Auditory';
roiNames{3} = 'SomatoSensory';
roiNames{4} = 'Retrosplenial';
roiNames{5} = 'Motor2';

% initialise matrices:
% for comparison of all trials that had a stimulus
anyC_percGO     = zeros(length(AL),nPerc,nROIs)*NaN;
anyC_percNOGO   = zeros(length(AL),nPerc,nROIs)*NaN;
% for comparison of all trials that had a stimulus, excluding the ones with equal contrasts
anyC_percGO_noEC   = zeros(length(AL),nPerc,nROIs)*NaN;
anyC_percNOGO_noEC = zeros(length(AL),nPerc,nROIs)*NaN;
anyC_percCorr_noEC = zeros(length(AL),nPerc,nROIs)*NaN;
anyC_percInc_noEC  = zeros(length(AL),nPerc,nROIs)*NaN;
% for high contrast trials (with contrasts on both sides as well)
highC_percGO    = zeros(length(AL),nPerc,nROIs)*NaN;
highC_percNOGO  = zeros(length(AL),nPerc,nROIs)*NaN;
highC_percCorr  = zeros(length(AL),nPerc,nROIs)*NaN;
highC_percInc   = zeros(length(AL),nPerc,nROIs)*NaN;
% for high contrast trials with contrast on one side only
highC_percGO_noCC    = zeros(length(AL),nPerc,nROIs)*NaN;
highC_percNOGO_noCC  = zeros(length(AL),nPerc,nROIs)*NaN;
highC_percCorr_noCC  = zeros(length(AL),nPerc,nROIs)*NaN;
highC_percInc_noCC   = zeros(length(AL),nPerc,nROIs)*NaN;
% for low contrast trials (with contrasts on both sides as well)
lowC_percGO    = zeros(length(AL),nPerc,nROIs)*NaN;
lowC_percNOGO  = zeros(length(AL),nPerc,nROIs)*NaN;
lowC_percCorr  = zeros(length(AL),nPerc,nROIs)*NaN;
lowC_percInc   = zeros(length(AL),nPerc,nROIs)*NaN;
% for high contrast trials with contrast on one side only
lowC_percGO_noCC    = zeros(length(AL),nPerc,nROIs)*NaN;
lowC_percNOGO_noCC  = zeros(length(AL),nPerc,nROIs)*NaN;
lowC_percCorr_noCC  = zeros(length(AL),nPerc,nROIs)*NaN;
lowC_percInc_noCC   = zeros(length(AL),nPerc,nROIs)*NaN;
% for zero contrast trials
zeroC_percCNG   = zeros(length(AL),nPerc,nROIs)*NaN;
zeroC_percFA    = zeros(length(AL),nPerc,nROIs)*NaN;


%% loop through experiments
for iie = 1:length(AL)
    
    Exps.animal     = AL{iie}.mousename;
    Exps.iseries    = AL{iie}.series;
    Exps.iexp       = AL{iie}.exp;
    
    expRef = strcat(Exps.iseries(1:4),'-',Exps.iseries(5:6),'-',Exps.iseries(7:8),...
        '_',Exps.iexp,'_',Exps.animal);
    
    if length(Exps.iexp) < 2        % only one experiment to analyse
        [b] = generateGenBlock(expRef, Exps);
        [U, V, t, SVDinfo]  = loadSVDfiles(b, true, true);
        [baseline, doPupil] = get_baseline(b, t, false);
        if baseline.firstTrialExcluded
            b.excludeFirstTrial = true;
        end
        ntr = b.completedTrials;
        if b.excludeFirstTrial
            ntr = ntr-1;
        end
    else
        [allExps,~,~,~] = concatenateExps(Exps.animal,expRef(1:10),str2num(Exps.iexp),true,true);
        b = allExps.block;
        U = allExps.U;
        V = allExps.V;
        t = allExps.t;
        if isequal(allExps.SVDinfo{1},allExps.SVDinfo{:})
            SVDinfo = allExps.SVDinfo{1};
        else        % this should be quite unlikely but just in case
            error('The separate experiments were processed differently and therefore should not be combined');
        end
        
        [baseline, doPupil] = get_baseline_concBlocks(b, t, doPupil);
        ntr = b.completedTrials;
        
        if isfield(b,'exps')      % for concatenated blocks
            b.iexp = num2str(b.exps);
        end
    end
    
    Fs = round(median(1/median(diff(t))));
    [baseline] = find_noMvmtBaseline(baseline, ntr, Fs);
    
    %% get behavioural variables
    [rGO, rNOGO, rCorr, rInc, rCNG, rGOc, rIncZC] = ...
        getChoiceNeglectInds(b,ntr);
    
    switch stimName
        case 'AudioVisual'
            if isfield(b,'excludeFirstTrial')
                if b.excludeFirstTrial
                    stim = b.stimuli.visContrasts(2:ntr+1,:);
                else
                    stim  = b.stimuli.visContrasts(1:ntr,:);
                end
            else
                stim  = b.stimuli.visContrasts(1:ntr,:);
            end
        otherwise
            if isfield(b,'excludeFirstTrial')
                if b.excludeFirstTrial
                    stim = b.stimuli(2:ntr+1,:);
                else
                    stim  = b.stimuli(1:ntr,:);
                end
            else
                stim  = b.stimuli(1:ntr,:);
            end
    end
    
    x1 = ismember(stim,[0 0],'rows');
    zeroCtrials = find(x1);
    clear x1
    
    equalCtrials    = find(diff(stim')==0);
    equalCNZtrials  = setdiff(equalCtrials,zeroCtrials);
    
    bothCtrials     = intersect(find(stim(:,1)>0),find(stim(:,2)>0));
    
    switch stimName
        case 'Auditory'
            highCond = 0.05;
        otherwise
            highCond = 0.5;
    end
    
    [highCtrials,~] = find(stim>=highCond); highCtrials = sort(highCtrials);
    highCtrialsNE   = setdiff(highCtrials,equalCNZtrials);
    highCtrialsNCC  = setdiff(highCtrials,bothCtrials);
    
    [lowCtrials_any,~]  = find(stim<highCond & stim~=0);
    lowCtrials      = setdiff(lowCtrials_any,highCtrials);
    lowCtrialsNE    = setdiff(lowCtrials,equalCNZtrials);
    lowCtrialsNCC   = setdiff(lowCtrials,bothCtrials);
    
    anyCtrials      = setdiff(1:size(stim,1),zeroCtrials);
    anyCtrialsNOEC  = setdiff(1:size(stim,1),equalCtrials);
    
    
    if ~isempty(intersect(zeroCtrials,lowCtrials)) && ~isempty(intersect(zeroCtrials,highCtrials)) ...
            && ~isempty(intersect(zeroCtrials,anyCtrials))
        error('something is wrong with the contrast indexing');
    end
    
    %% load or compute powermap
    if exist(fullfile(myDirectoryWithTheData,'BaselinePower', ...
            strcat(b.animal,'_',b.iseries,'_',b.iexp,'_',...
            num2str(frb(1)),'to',num2str(frb(2)),'Hz','_PowerMap_noMvmt.mat')));
        
        disp('loading power map...');
        load(fullfile(myDirectoryWithTheData,'BaselinePower', ...
            strcat(b.animal,'_',b.iseries,'_',b.iexp,'_',...
            num2str(frb(1)),'to',num2str(frb(2)),'Hz','_PowerMap_noMvmt.mat')),'Pmap');
        disp('done.');
        
    else
        
        [U, V, t, SVDinfo]  = loadSVDfiles(b, true, true);
        [Pmap] = computePowerMap_noMvmt(b, ntr, U, V, t, frb, baseline, shortbl);
        
    end
    
    if size(Pmap,3) < ntr
        ntr = size(Pmap,3);
    end
    
    %%
    
    load(fullfile(myDirectoryWithTheData,'ROIPixelSelections',...
        strcat(b.animal,'_',b.iseries,'_',b.iexp,'_PixelPerROI.mat')));
    
    for iir = 1:nROIs
        if isfinite(pix(1,iir))
            thisPixP = log10(squeeze(Pmap(pix(1,iir), pix(2,iir),:)));
            powerPercentiles = prctile(thisPixP,percentiles);
            
            for prc = 1:nPerc
                thisPrcTrials = find(thisPixP>=powerPercentiles(prc) & ...
                    thisPixP<powerPercentiles(prc+1));
                
                for nc = 1:7    % loop through the stimulus conditions
                    switch nc
                        case 1
                            theseTrials = intersect(thisPrcTrials,anyCtrials);
                            if ~isempty(theseTrials)
                                if length(intersect(theseTrials,rGO))+length(intersect(theseTrials,rNOGO)) ...
                                        ~= length(theseTrials)
                                    error('something does not add up');
                                end
                                anyC_percGO(iie,prc,iir) = ...
                                    length(intersect(theseTrials,rGO))/length(theseTrials);
                                anyC_percNOGO(iie,prc,iir) = ...
                                    length(intersect(theseTrials,rNOGO))/length(theseTrials);
                            end
                            
                        case 3
                            theseTrials = intersect(thisPrcTrials,highCtrialsNE);
                            if ~isempty(theseTrials)
                                if length(intersect(theseTrials,rCorr))+length(intersect(theseTrials,rInc)) ...
                                        + length(intersect(theseTrials,rNOGO)) ~= length(theseTrials)
                                    error('something does not add up');
                                end
                                highC_percGO(iie,prc,iir) = ...
                                    length(intersect(theseTrials,rCorr))/length(theseTrials);
                                highC_percNOGO(iie,prc,iir) = ...
                                    length(intersect(theseTrials,rNOGO))/length(theseTrials);
                                
                                theseTrialsCorr = intersect(theseTrials,rCorr);
                                theseTrialsInc  = intersect(theseTrials,rInc);
                                theseTrialsSum  = length(theseTrialsCorr)+length(theseTrialsInc);
                                if theseTrialsSum > 0
                                    highC_percCorr(iie,prc,iir) = ...
                                        length(theseTrialsCorr)/theseTrialsSum;
                                    highC_percInc(iie,prc,iir) = ...
                                        length(theseTrialsInc)/theseTrialsSum;
                                end
                            end
                            
                        case 5
                            theseTrials = intersect(thisPrcTrials,lowCtrialsNE);
                            if ~isempty(theseTrials)
                                if length(intersect(theseTrials,rCorr))+length(intersect(theseTrials,rInc)) ...
                                        + length(intersect(theseTrials,rNOGO)) ~= length(theseTrials)
                                    error('something does not add up');
                                end
                                lowC_percGO(iie,prc,iir) = ...
                                    length(intersect(theseTrials,rCorr))/length(theseTrials);
                                lowC_percNOGO(iie,prc,iir) = ...
                                    length(intersect(theseTrials,rNOGO))/length(theseTrials);
                                
                                theseTrialsCorr = intersect(theseTrials,rCorr);
                                theseTrialsInc  = intersect(theseTrials,rInc);
                                theseTrialsSum  = length(theseTrialsCorr)+length(theseTrialsInc);
                                if theseTrialsSum > 0
                                    lowC_percCorr(iie,prc,iir) = ...
                                        length(theseTrialsCorr)/theseTrialsSum;
                                    lowC_percInc(iie,prc,iir) = ...
                                        length(theseTrialsInc)/theseTrialsSum;
                                end
                            end
                            
                        case 7
                            theseTrials = intersect(thisPrcTrials,zeroCtrials);
                            if ~isempty(theseTrials)
                                if length(intersect(theseTrials,rCNG))+length(intersect(theseTrials,rIncZC)) ...
                                        ~= length(theseTrials)
                                    error('something does not add up');
                                end
                                zeroC_percCNG(iie,prc,iir) = ...
                                    length(intersect(theseTrials,rCNG))/length(theseTrials);
                                zeroC_percFA(iie,prc,iir) = ...
                                    length(intersect(theseTrials,rIncZC))/length(theseTrials);
                            end
                            
                        case 2
                            theseTrials = intersect(thisPrcTrials,anyCtrialsNOEC);
                            if ~isempty(theseTrials)
                                if length(intersect(theseTrials,rCorr)) + length(intersect(theseTrials,rInc)) ...
                                        + length(intersect(theseTrials,rNOGO)) ~= length(theseTrials)
                                    error('something does not add up');
                                end
                                anyC_percGO_noEC(iie,prc,iir)   = ...
                                    length(intersect(theseTrials,rGO))/length(theseTrials);
                                anyC_percNOGO_noEC(iie,prc,iir)   = ...
                                    length(intersect(theseTrials,rNOGO))/length(theseTrials);
                                
                                theseTrialsCorr = intersect(theseTrials,rCorr);
                                theseTrialsInc  = intersect(theseTrials,rInc);
                                theseTrialsSum  = length(theseTrialsCorr)+length(theseTrialsInc);
                                if theseTrialsSum > 0
                                    anyC_percCorr_noEC(iie,prc,iir) = ...
                                        length(theseTrialsCorr)/theseTrialsSum;
                                    anyC_percInc_noEC(iie,prc,iir) = ...
                                        length(theseTrialsInc)/theseTrialsSum;
                                end
                            end
                            
                        case 4
                            theseTrials = intersect(thisPrcTrials,highCtrialsNCC);
                            if ~isempty(theseTrials)
                                if length(intersect(theseTrials,rCorr))+length(intersect(theseTrials,rInc)) ...
                                        + length(intersect(theseTrials,rNOGO)) ~= length(theseTrials)
                                    error('something does not add up');
                                end
                                highC_percGO_noCC(iie,prc,iir) = ...
                                    length(intersect(theseTrials,rCorr))/length(theseTrials);
                                highC_percNOGO_noCC(iie,prc,iir) = ...
                                    length(intersect(theseTrials,rNOGO))/length(theseTrials);
                                
                                theseTrialsCorr = intersect(theseTrials,rCorr);
                                theseTrialsInc  = intersect(theseTrials,rInc);
                                theseTrialsSum  = length(theseTrialsCorr)+length(theseTrialsInc);
                                if theseTrialsSum
                                    highC_percCorr_noCC(iie,prc,iir) = ...
                                        length(theseTrialsCorr)/theseTrialsSum;
                                    highC_percInc_noCC(iie,prc,iir) = ...
                                        length(theseTrialsInc)/theseTrialsSum;
                                end
                            end
                            
                        case 6
                            theseTrials = intersect(thisPrcTrials,lowCtrialsNCC);
                            if ~isempty(theseTrials)
                                if length(intersect(theseTrials,rCorr))+length(intersect(theseTrials,rInc)) ...
                                        + length(intersect(theseTrials,rNOGO)) ~= length(theseTrials)
                                    error('something does not add up');
                                end
                                lowC_percGO_noCC(iie,prc,iir) = ...
                                    length(intersect(theseTrials,rCorr))/length(theseTrials);
                                lowC_percNOGO_noCC(iie,prc,iir) = ...
                                    length(intersect(theseTrials,rNOGO))/length(theseTrials);
                                
                                theseTrialsCorr = intersect(theseTrials,rCorr);
                                theseTrialsInc  = intersect(theseTrials,rInc);
                                theseTrialsSum  = length(theseTrialsCorr)+length(theseTrialsInc);
                                if theseTrialsSum > 0
                                    lowC_percCorr_noCC(iie,prc,iir) = ...
                                        length(theseTrialsCorr)/theseTrialsSum;
                                    lowC_percInc_noCC(iie,prc,iir) = ...
                                        length(theseTrialsInc)/theseTrialsSum;
                                end
                            end
                            
                    end
                end
            end
            
        end
    end
    
    
end


%% save
saveDir = myDirectoryWithTheData;
fname = strcat(stimName,'_summary_behaviourPerPower_',num2str(frb(1)),'to',num2str(frb(2)),'Hz',...
    '_with',num2str(nPerc),'percentiles');
save(fullfile(saveDir,fname),...
    'anyC_percGO','anyC_percNOGO',...
    'anyC_percGO_noEC','anyC_percNOGO_noEC','anyC_percCorr_noEC','anyC_percInc_noEC',...
    'highC_percGO','highC_percNOGO','highC_percCorr','highC_percInc',...
    'highC_percGO_noCC','highC_percNOGO_noCC','highC_percCorr_noCC','highC_percInc_noCC',...
    'lowC_percGO','lowC_percNOGO','lowC_percCorr','lowC_percInc',...
    'lowC_percGO_noCC','lowC_percNOGO_noCC','lowC_percCorr_noCC','lowC_percInc_noCC',...
    'zeroC_percCNG','zeroC_percFA',...
    'percentiles','nPerc','nROIs','roiNames','frb');


%% plot

whichMice = 'all'; % all, EJ (unilateral), NS (bilateral), or indicator
switch whichMice
    case 'all'
        expInds = [1:58];
    case 'EJ'
        expInds = [1:37];
    case 'NS'
        expInds = [38:58];
    case 'GCamp6s'
        expInds = [1,20:22,41:56];
    case 'GCamp6f'
        expInds = [2:19,23:40,57:58];
end

figure;
for cc = 1:6
    for rr = 1:nROIs
        % plot the different ROIs in different columns
        switch rr
            case 1
                x = 0;
            case 2
                x = 6;
            case 3
                x = 12;
            case 4
                x = 18;
            case 5
                x = 24;
        end
        subplot(5,6,cc+x)
        switch cc
            case 1
                shadedErrorBar(percentiles(2:end),nanmean(anyC_percGO(expInds,:,rr)),...
                    nansem(anyC_percGO(expInds,:,rr),[],1),'-y',1); hold on;
                shadedErrorBar(percentiles(2:end),nanmean(anyC_percNOGO(expInds,:,rr)),...
                    nansem(anyC_percNOGO(expInds,:,rr),[],1),'-k',1);
                title(['Any stim']);
            case 2
                shadedErrorBar(percentiles(2:end),nanmean(highC_percGO(expInds,:,rr)),...
                    nansem(highC_percGO(expInds,:,rr),[],1),'--y',1); hold on;
                shadedErrorBar(percentiles(2:end),nanmean(highC_percNOGO(expInds,:,rr)),...
                    nansem(highC_percNOGO(expInds,:,rr),[],1),'--k',1);
                shadedErrorBar(percentiles(2:end),nanmean(lowC_percGO(expInds,:,rr)),...
                    nansem(lowC_percGO(expInds,:,rr),[],1),'-.y',1); hold on;
                shadedErrorBar(percentiles(2:end),nanmean(lowC_percNOGO(expInds,:,rr)),...
                    nansem(lowC_percNOGO(expInds,:,rr),[],1),'-.k',1);
                title(['High + low (all)']);
                
                
            case 3
                shadedErrorBar(percentiles(2:end),nanmean(highC_percGO_noCC(expInds,:,rr)),...
                    nansem(highC_percGO_noCC(expInds,:,rr),[],1),'--y',1); hold on;
                shadedErrorBar(percentiles(2:end),nanmean(highC_percNOGO_noCC(expInds,:,rr)),...
                    nansem(highC_percNOGO_noCC(expInds,:,rr),[],1),'--k',1);
                shadedErrorBar(percentiles(2:end),nanmean(lowC_percGO_noCC(expInds,:,rr)),...
                    nansem(lowC_percGO_noCC(expInds,:,rr),[],1),'-.y',1); hold on;
                shadedErrorBar(percentiles(2:end),nanmean(lowC_percNOGO_noCC(expInds,:,rr)),...
                    nansem(lowC_percNOGO_noCC(expInds,:,rr),[],1),'-.k',1);
                title(['High + low (single C)']);
                
            case 4
                shadedErrorBar(percentiles(2:end),nanmean(anyC_percCorr_noEC(expInds,:,rr)),...
                    nansem(anyC_percCorr_noEC(expInds,:,rr),[],1),'-g',1); hold on;
                shadedErrorBar(percentiles(2:end),nanmean(anyC_percInc_noEC(expInds,:,rr)),...
                    nansem(anyC_percInc_noEC(expInds,:,rr),[],1),'-r',1);
                
                shadedErrorBar(percentiles(2:end),nanmean(highC_percCorr(expInds,:,rr)),...
                    nansem(highC_percCorr(expInds,:,rr),[],1),'--g',1); hold on;
                shadedErrorBar(percentiles(2:end),nanmean(highC_percInc(expInds,:,rr)),...
                    nansem(highC_percInc(expInds,:,rr),[],1),'--r',1);
                
                shadedErrorBar(percentiles(2:end),nanmean(lowC_percCorr(expInds,:,rr)),...
                    nansem(lowC_percCorr(expInds,:,rr),[],1),'-.g',1); hold on;
                shadedErrorBar(percentiles(2:end),nanmean(lowC_percInc(expInds,:,rr)),...
                    nansem(lowC_percInc(expInds,:,rr),[],1),'-.r',1);
                title(['Any stim (noEC) + high + low']);
                
            case 5
                shadedErrorBar(percentiles(2:end),nanmean(highC_percCorr_noCC(expInds,:,rr)),...
                    nansem(highC_percCorr_noCC(expInds,:,rr),[],1),'--g',1); hold on;
                shadedErrorBar(percentiles(2:end),nanmean(highC_percInc_noCC(expInds,:,rr)),...
                    nansem(highC_percInc_noCC(expInds,:,rr),[],1),'--r',1);
                
                shadedErrorBar(percentiles(2:end),nanmean(lowC_percCorr_noCC(expInds,:,rr)),...
                    nansem(lowC_percCorr_noCC(expInds,:,rr),[],1),'-.g',1); hold on;
                shadedErrorBar(percentiles(2:end),nanmean(lowC_percInc_noCC(expInds,:,rr)),...
                    nansem(lowC_percInc_noCC(expInds,:,rr),[],1),'-.r',1);
                title(['High + low (single C)']);
                
            case 6
                shadedErrorBar(percentiles(2:end),nanmean(zeroC_percCNG(expInds,:,rr)),...
                    nansem(zeroC_percCNG(expInds,:,rr),[],1),'-c',1); hold on;
                shadedErrorBar(percentiles(2:end),nanmean(zeroC_percFA(expInds,:,rr)),...
                    nansem(zeroC_percFA(expInds,:,rr),[],1),'-m',1);
                title(['Zero contrast, ' roiNames{rr}]);
        end
        
        xlabel([num2str(frb(1)) '-' num2str(frb(2)) 'Hz P percentiles']);
        ylabel('Resp. (%)');
        ylim([0 1]);
        axis square
        set(gca, 'FontSize', 12)
    end
    
end


%% figure to compare each area

figure;
for sp = 1:12
    switch sp
        case 1
            thisPerc = anyC_percGO;
            ttl = 'All C Choice';
        case 2
            thisPerc = highC_percGO;
            ttl = 'High C Choice';
        case 3
            thisPerc = lowC_percGO;
            ttl = 'Low C Choice';
        case 4
            thisPerc = highC_percCorr;
            ttl = 'High C Correct';
        case 5
            thisPerc = lowC_percCorr;
            ttl = 'Low C Correct';
        case 6
            thisPerc = zeroC_percCNG;
            ttl = 'Zero C Corr Reject';
        case 7
            thisPerc = anyC_percNOGO;
            ttl = 'All C Neglect';
        case 8
            thisPerc = highC_percNOGO;
            ttl = 'High C Neglect';
        case 9
            thisPerc = lowC_percNOGO;
            ttl = 'Low C Neglect';
        case 10
            thisPerc = highC_percInc;
            ttl = 'High C Incorrect';
        case 11
            thisPerc = lowC_percInc;
            ttl = 'Low C Incorrect';
        case 12
            thisPerc = zeroC_percFA;
            ttl = 'Zero C False Alarm';
    end
    
    subplot(2,6,sp)
    for rr = [1 3] %1:5
        switch rr
            case 1
                col = 'b';
            case 2
                col = 'r';
            case 3
                col = 'g';
            case 4
                col = 'k';
            case 5
                col = 'y';
        end
        shadedErrorBar(percentiles(2:end),nanmean(thisPerc(expInds,:,rr)),...
            nansem(thisPerc(expInds,:,rr),[],1),strcat('-',col),1); hold on;
        %         plot(percentiles(2:end),nanmean(thisPerc(expInds,:,rr)),col); hold on
        title(ttl);
        xlabel([num2str(frb(1)) '-' num2str(frb(2)) 'Hz P percentiles']);
        ylabel('Resp. (%)');
        ylim([0 1]);
        axis square
        set(gca, 'FontSize', 12)
    end
end

