% script with the code that generated the visual physiology during behaviour example panels
% in Figure 2
%
% written by Elina Jacobs, UCL Cortexlab

% set directories here
thisDir = myDirectoryWithTheData;

useDFF = false;

Exps.animal     = 'Cori';
Exps.iseries    = '20161208';
Exps.iexp       = '1';

expRef = strcat(Exps.iseries(1:4),'-',Exps.iseries(5:6),'-',Exps.iseries(7:8),...
    '_',Exps.iexp,'_',Exps.animal);

[b] = generateGenBlock(expRef, Exps);
[U, V, t, SVDinfo]  = loadSVDfiles(b, true, true);

if useDFF
    [dffU, dffV] = applyDFF(U,V,b);
end

[baseline, ~] = get_baseline(b, t, false);
shortbl = find(diff(baseline.OnOffTimes')<1);

ntr = b.completedTrials;
if b.excludeFirstTrial
    ntr = ntr-1;
end
tEt   = [b.evts.endTrialTimes(1:ntr)]./60;
tSt   = [b.evts.newTrialTimes(1:ntr)]./60;

Fs  = round(1/median(diff(t)));

frb = [3 6];

nSV = size(U,3);

%% generate stimulus triggered average image

% load ROI pixel coordinates
load(fullfile(myDirectoryWithTheData,'ROIPixelSelections',...
    strcat(b.animal,'_',b.iseries,'_',b.iexp,'_PixelPerROI.mat')));

if useDFF
    Vr = V;         % save non-dff V under different name so can call it later again for spectral analysis
    V  = dffV;
    Ur = U;         % save non-dff U under different name
    U  = dffU;
end

tS = [0.05 0.07];   

stimInds = find(b.stimuli(:,2)>0.55);
while stimInds(end) > ntr
    stimInds(end) = [];
end

stimTimes = b.evts.stimuliOnTimes(1:ntr);
% % calculate baseline at time 0 (stimulus onset)
baselineV = zeros(size(V,1), ntr, 1);
for s = 1:nSV
    baselineV(s,:) = interp1(t, V(s,:), stimTimes);
end

% calculate stimulus response in time window specified by tS
winSampsStim = tS(1):(1/Fs):tS(2);
stimInt = bsxfun(@plus, stimTimes', winSampsStim);

stimV = zeros(size(V,1), ntr, length(winSampsStim));
for s = 1:nSV
    stimV(s,:,:) = interp1(t, V(s,:), stimInt);
end
stimV       = bsxfun(@minus, stimV, baselineV);
stimAvgV    = nanmean(stimV,3);
stimRespAvg = reshape(reshape(U,size(U,1)*size(U,2),size(U,3))*stimAvgV,size(U,1),size(U,2),ntr);

SRstim = mean(stimRespAvg(:,:,stimInds),3);
SRall  = mean(stimRespAvg,3);
SR = SRstim - SRall;
% SR = SRstim;
if useDFF
    xl = max(abs(SR(:)))*0.1;           % rounding loeads to xl = 0 as values from dff so much smaller
else
    xl = round(max(abs(SR(:)))*0.9);
end

% load powermap to get mask
disp('loading power map...');
load(fullfile(EJDirs.Data.Im,Exps.animal,'BaselinePower', ...
    strcat(b.animal,'_',b.iseries,'_',b.iexp,'_',...
    num2str(frb(1)),'to',num2str(frb(2)),'Hz','_PowerMap.mat')),'Pmap');
disp('done.');

Pmap(:,:,shortbl) = NaN;

meanP = log10(nanmean(Pmap,3));

%
cm = RedWhiteBlue; cm = flipud(cm);

cAxLims = [-xl xl; ...
    round(prctile(meanP(:),25),2) round(prctile(meanP(:),75),2)];

stimRespMap = createRGBcomposite(SR,meanP,cAxLims);

figure;
image(stimRespMap);
hold on; axis off; axis image;
plot(pix(2,1),pix(1,1),'k.','MarkerSize',32);
if useDFF
    title(['Stimulus response at ' num2str(tS(1)*1000) '-' num2str(tS(2)*1000) ' ms (from dff)']);
else
    title(['Stimulus response at ' num2str(tS(1)*1000) '-' num2str(tS(2)*1000) ' ms']);
end

% how plot in Figure 2 was created:
% compute SR from 0.05-0.07s, then subtract SR from 0.07-0.08s, subtract
% the first from the latter -> that is the stimulus response
% (and SR is SRstim - SRall, in each case)

%% example traces

[baseline] = find_noMvmtBaseline(baseline, ntr, Fs);
% get visual cortex trace
% visual cortex ROI
load(fullfile(myDirectoryWithTheData,'ROIPixelSelections',...
    strcat(b.animal,'_',b.iseries,'_',b.iexp,'_PixelPerROI.mat')));
pixU = squeeze(U(pix(1,1),pix(2,1),:))';
pix_trace = pixU*V;

wheel_t = interp1(b.wheel.Times, b.wheel.Values, t);

NOGOind = 397; % 397
GOind   = 218; % 218

figure;
suptitle([b.animal ' ' b.iseries ' ' b.iexp]);
for cc = 1:2
    switch cc
        case 1
            thisInd = NOGOind;
            tlt = 'Miss';
        case 2
            thisInd = GOind;
            tlt = 'Choice';
    end
    subplot(2,2,cc)
    cueOnset       = b.evts.cueOnTimes(2,thisInd);
    thisStimOnset  = b.evts.stimuliOnTimes(thisInd);
    thisStimOffset = b.evts.stimuliOffTimes(thisInd);
    noMvmtOnset    = t(baseline.noMvmtFrames(thisInd,1));
    noMvmtTime  = thisStimOnset - noMvmtOnset;
    
    rectangle('Position',[noMvmtOnset -2500 noMvmtTime 5000],...
        'FaceColor',[1 1 0.5],'EdgeColor','w'); hold on;
    plot([thisStimOnset thisStimOnset],...
        [-2500 2000],'Color',[0.8 0.8 0.8],'LineWidth',0.5);
    plot([cueOnset cueOnset],...
        [-2500 2000],'k:','LineWidth',2);
    plot([thisStimOffset thisStimOffset],...
        [-2500 2000],'Color',[0.8 0.8 0.8],'LineWidth',0.5);
    thisRV = b.evts.responseValues(thisInd);
    thisFT = b.evts.feedbackTimes(thisInd);
    % plot a dot to indicate trial outcome
    if thisRV == 0
        plot([thisFT thisFT],[-2500 2000],'--','Color',[0.3 0.3 0.3],'LineWidth',3);
    else
        thisFV = b.evts.feedbackValues(thisInd);
        if thisFV > 0
            plot([thisFT thisFT],[-2500 2000],'--','Color',[0 0.9 0.4],'LineWidth',3);
        else
            plot([thisFT thisFT],[-2500 2000],'--','Color',[0.64 0.08 0.18],'LineWidth',3);
        end
    end
    plot(t,pix_trace,'b','LineWidth',2);
    %     plot(t,bp_pixTrace,'c:');
    ylim([-2500 1800]);
    ylabel('F');
    xlim([cueOnset-2.5, cueOnset+2.5]);
    xlabel('Time (s)');
    title(tlt);
    set(gca, 'FontSize', 24)
    set(get(gca,'YLabel'),'Rotation',0)
    
    dist = abs(b.wheel.Times - cueOnset+1);
    idx = find(dist == min(dist));
    wheelVal = b.wheel.Values(idx);
    
    % plot wheel trace
    subplot(2,2,cc+2)
    rectangle('Position',[noMvmtOnset 12350 noMvmtTime 5000],...
        'FaceColor',[1 1 0.5],'EdgeColor','w'); hold on;
    plot([thisStimOnset thisStimOnset],...
        [12350 18000],'Color',[0.8 0.8 0.8],'LineWidth',0.5);
    plot([cueOnset cueOnset],...
        [12350 18000],'k:','LineWidth',2);
    plot([thisStimOffset thisStimOffset],...
        [12350 18000],'Color',[0.8 0.8 0.8],'LineWidth',0.5);
    thisRV = b.evts.responseValues(thisInd);
    thisFT = b.evts.feedbackTimes(thisInd);
    % plot a dot to indicate trial outcome
    if thisRV == 0
        plot([thisFT thisFT],[12350 18000],'--','Color',[0.3 0.3 0.3],'LineWidth',3);
    else
        thisFV = b.evts.feedbackValues(thisInd);
        if thisFV > 0
            plot([thisFT thisFT],[12350 18000],'--','Color',[0 0.9 0.4],'LineWidth',3);
        else
            plot([thisFT thisFT],[12350 18000],'--','Color',[0.64 0.08 0.18],'LineWidth',3);
        end
    end
    
    plot(b.wheel.Times,b.wheel.Values,'Color',[0.85 0.6 0.4],'LineWidth',2);
    ylim([wheelVal-200 wheelVal+800]);
    ylabel('wheel trace (au)');
    xlim([cueOnset-2.5, cueOnset+2.5]);
    xlabel('Time (s)');
    title(tlt);
    set(gca, 'FontSize', 24);
    set(gca, 'YTickLabel',[]);
    set(get(gca,'YLabel'),'Rotation',0)
    
end