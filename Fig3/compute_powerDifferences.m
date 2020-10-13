% script that calculates power differences during quiescent periods across 
% experiments specified by the expList

% written by Elina Jacobs, UCL Cortexlab

% set directories here
thisDir = myDirectoryWithTheData;

stimName = 'VisualALL';

load(fullfile(myDirectoryWithTheData,...
    strcat(stimName,'_','expList_forPowerAnalysis.mat')));

AL = expList_forPowerAnalysis; clear expList_forPowerAnalysis

%% initialisations
frb         = [3,6];    % frequency band of interest

% say whether power to be calculated in ROI or as average across window
pr          = 'ROIs';   
switch pr
    case 'average'
        nROIs       = 1;
        roiNames{1} = 'windowAverage';
        
    case 'ROIs'
        % edited on 2017/10/17 to include more ROIs
        switch stimName
            case 'Visual'
                nROIs = 5;
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
        
end

minBL       = 0.7;      % minimum baseline length in seconds
minTrials   = 10;       % minimum number of trials per behavioural condition 

whichCond   = 'Choice_Miss';
% Choice_Miss
% Correct_Miss
% Incorrect_Miss
% Correct_Incorrect
% FalseAlarm_CorrectReject
% CorrectReject_Miss
% Correct_Incorrect_Reward
% CorrectReject_Incorrect_Reward

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
    
    %     shortbl = find(diff(baseline.OnOffTimes')<0.7);
    Fs = round(median(1/median(diff(t))));
    [baseline] = find_noMvmtBaseline(baseline, ntr, Fs);
    
    shortbl = find(diff(baseline.noMvmtFrames')<minBL*Fs);
    disp(Exps);
    disp(strcat(num2str(round(length(shortbl)/ntr*100)),'% of trials have short baseline'));

    
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
        
    
    %% get behavioural variables
    
    % if we're interesed in the reward, ie post trial and post movement
    % power difference, exclude the last trial
    switch whichCond
        case 'Correct_Incorrect_Reward'
            ntr = ntr-1;
        case 'CorrectReject_Incorrect_Reward'
            ntr = ntr-1;
    end
    
    [rGO, rNOGO, rGoCorr, rGoInc, rCNG, rGOc, rInZC] = getChoiceNeglectInds(b,ntr,stimName);
    % rGO includes false alarms; so this is all trials where the animals responded regardless of whether there was a contrast or not
    % rGOc specifies the trials where a go response was made when there was a contrast
    % rInZC specifies the trials where a go response was made when there was no contrast
    
    switch stimName
        case 'Auditory'
            switch AL{iie}.mousename
                case 'EJ015'
                    rT = b.evts.responseTimes - b.evts.interactiveOnT;
                    zz = find(rT>7.5);
                    rNOGO=zz;
                    rGO(zz) = [];
            end
    end
    
    switch whichCond
        case 'Choice_Miss'
            cond1 = rGOc;           
            cond2 = rNOGO;
        case 'Correct_Incorrect'
            cond1 = rGoCorr;
            cond2 = rGoInc;
        case 'FalseAlarm_CorrectReject';
            cond1 = rInZC;
            cond2 = rCNG;
        case 'CorrectReject_Miss';
            cond1 = rCNG;
            cond2 = rNOGO;
        case 'Correct_Miss'
            cond1 = rGoCorr;
            cond2 = rNOGO;
        case 'Incorrect_Miss'
            cond1 = rGoInc;
            cond2 = rNOGO;
        case 'Correct_Incorrect_Reward'
            cond1 = rGoCorr+1;
            cond2 = rGoInc+1;
        case 'CorrectReject_Incorrect_Reward'
            cond1 = rCNG+1;
            cond2 = rGoInc+1;
    end
            
    
    %% pick pixels for power analysis
    
    switch pr
        case 'ROIs'
            % load pixel coordinates
            load(fullfile(myDirectoryWithTheData,'ROIPixelSelections',...
                strcat(b.animal,'_',b.iseries,'_',b.iexp,'_PixelPerROI.mat')));
            
            for iir = 1:nROIs
                
                if isfinite(pix(1,iir))
                    %                     thisPixP = squeeze(log10(Pmap(pix(1,iir), pix(2,iir),:)));
                    thisPixP = squeeze(Pmap(pix(1,iir), pix(2,iir),:));
                    
                    if length(cond1) >= minTrials && length(cond2) >= minTrials
                        PPP(1,iir,iie)  =  log10(nanmean(thisPixP(cond1)));
                        SSEM(1,iir,iie) =  log10(nansem(thisPixP(cond1),[],1));
                        
                        PPP(2,iir,iie)  =  log10(nanmean(thisPixP(cond2)));
                        SSEM(2,iir,iie) =  log10(nansem(thisPixP(cond2),[],1));
                        
                        pvals(1,iir,iie) = ranksum(thisPixP(cond1),thisPixP(cond2));
                    else
                        PPP(1:2,iir,iie)  = NaN;
                        SSEM(1:2,iir,iie) = NaN;
                        pvals(1,iir,iie)  = NaN;
                    end
                    
                else
                    PPP(1:2,iir,iie)  = NaN;
                    SSEM(1:2,iir,iie) = NaN;
                    pvals(1,iir,iie)  = NaN;
                end
            end
            
            save(fullfile(myDirectoryWithTheData,...
                strcat(stimName,'_',num2str(frb(1)),'to',num2str(frb(2)),'Hz_Power_perROI_',whichCond,'.mat')),'PPP','SSEM','pvals');
            
            close all
            
        case 'average'
            
            % this bit will set pixels outside the brain to NaN
            meanP = nanmean(Pmap,3);
            prc10 = find(Pmap<=prctile(meanP(:),20));
            maskPmap = Pmap;
            maskPmap(prc10) = NaN;
            
            flatPmap = reshape(maskPmap,size(maskPmap,1)*size(maskPmap,2),size(maskPmap,3));
            xx = find(flatPmap==0);
            flatPmap(xx) = NaN; clear xx
            avgP = nanmean(flatPmap);
            
            PPP(1,iie)  =  log10(nanmean(avgP(cond1)));
            SSEM(1,iie) =  log10(nansem(avgP(cond1),[],2));
            
            PPP(2,iie)  =  log10(nanmean(avgP(cond2)));
            SSEM(2,iie) =  log10(nansem(avgP(cond2),[],2));
            
            if ~isempty(cond1) && ~isempty(cond2)
                pvals(iie) = ranksum(avgP(cond1),avgP(cond2));
            else pvals(iie) = NaN;
            end
            
            save(fullfile(myDirectoryWithTheData,...
                strcat(stimName,'_Power_avg_',whichCond,'.mat')),'PPP','pvals');
    end
    
end

