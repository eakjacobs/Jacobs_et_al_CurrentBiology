% examplePowerDifferenceMaps
% script that computes the power difference maps
% written by Elina Jacobs, UCL Cortexlab

%% initialisations
% set directory with the data
thisDir = myDirectoryWithTheData;

stimName = 'VisualALL'; 
% specify the behavioural experiment type: Visual, Auditory,
% AudioVisual (refers to the auditory distractor task)
% VisualALL contains both 2AFC, 2AUC and Contrast Discrimination data,
% Visual2ANFC excludes the 2AFC tasks

whichMap = 'ChoiceMiss';    
% CorrectIncorrect
% FaCorrRej (False Alarm - Correct Reject)
% CorrRejMiss

% datasets used as examples in Figure 3
% unilateral imaging dataset:
Exps.animal      = 'EJ013';
Exps.iseries     = '20160601';
Exps.iexp        = '1';
% bilateral imaging dataset:
Exps.animal     = 'Cori';
Exps.iseries    = '20161208';
Exps.iexp       = '1';

% datasets used as examples in Figure 4
% unilateral imaging dataset:
Exps.animal   = 'EJ009';
Exps.iseries  = '20150729';
Exps.iexp     = '1';
% bilateral example is the same is in Figure 3
                
% datasets used as examples in Figure 5
% unilateral example is the same is in Figure 3
% bilateral example is the same is in Figure 3

% datasets used as examples in Figure 6
% unilateral example is the same is in Figure 4
% bilateral imaging dataset:
Exps.animal     = 'Cori';
Exps.iseries    = '20161210';
Exps.iexp       = '1';

% dataset used as example in Figure 7
Exps.animal   = 'EJ011';
Exps.iseries  = '20161004';
Exps.iexp     = '3';

% dataset used as example in Supplementary Figure 7
Exps.animal   = 'EJ007';
Exps.iseries  = '20151209';
Exps.iexp     = '1';

% comment out the irrelevant experiments

expRef = strcat(Exps.iseries(1:4),'-',Exps.iseries(5:6),'-',Exps.iseries(7:8),...
    '_',Exps.iexp,'_',Exps.animal);

% load behavioural data
[b] = generateGenBlock(expRef, Exps);

% if powermap hasn't been computed yet, run computePowerMap_noMvmt
% load powermap
frb         = [3,6];
disp('loading power map...');
load(fullfile(myDirectoryWithTheData,'BaselinePower', ...
    strcat(b.animal,'_',b.iseries,'_',b.iexp,'_',...
    num2str(frb(1)),'to',num2str(frb(2)),'Hz','_PowerMap_noMvmt.mat')),'Pmap');
disp('done.');

% load ROI pixel coordinates
load(fullfile(myDirectoryWithTheData,'ROIPixelSelections',...
    strcat(b.animal,'_',b.iseries,'_',b.iexp,'_PixelPerROI.mat')));

% specify some plot details
cm = RedWhiteBlue; cm = flipud(cm);
% changed on 20180127 from 25-75 prctile to 10-90
cAxLims = [-0.75 0.75; ...
    round(prctile(meanP(:),20),2) round(prctile(meanP(:),25),2)];

%% computations
% calculate mean power to use as mask
meanP = log10(nanmean(Pmap,3));

% get trial indices with the different behavioural responses
[rGO, rNOGO, rGoCorr, rGoInc, rCNG, rGOc, rInZC] = getChoiceNeglectInds_separateZC(b,ntr,stimName);

switch whichMap
    case 'ChoiceMiss'
        Pmap_Choice = nanmean(Pmap(:,:,rGO),3);
        Pmap_Miss = nanmean(Pmap(:,:,rNOGO),3);
        PdiffMap = log10(Pmap_Choice) - log10(Pmap_Miss);
    case 'CorrectIncorrect'
        Pmap_Correct = nanmean(Pmap(:,:,rGoCorr),3);
        Pmap_Incorrect = nanmean(Pmap(:,:,rGoInc),3);
        PdiffMap = log10(Pmap_Correct) - log10(Pmap_Incorrect);
    case 'FaCorrRej'
        Pmap_FA = nanmean(Pmap(:,:,rInZC),3);
        Pmap_CorrRej = nanmean(Pmap(:,:,rCNG),3);
        PdiffMap = log10(Pmap_FA) - log10(Pmap_CorrRej);
    case 'CorrRejMiss'
        Pmap_CorrRej = nanmean(Pmap(:,:,rCNG),3);
        Pmap_Miss = nanmean(Pmap(:,:,rNOGO),3);
        PdiffMap = log10(Pmap_CorrRej) - log10(Pmap_Miss);
        
end

map1 = createRGBcomposite(PdiffMap,meanP,cAxLims);

figure;
image(map1); axis off; axis image; hold on;
% plot ROI pixel points
for iip = 1:size(pix,2)
    plot(pix(2,iip),pix(1,iip),'k.','MarkerSize',18);
end
title(whichMap)