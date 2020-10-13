% script that generates powerspectrum ratio in visual cortex for every
% single dataset; then plots average with SEM of that
%
% written by Elina Jacobs, UCL Cortexlab

% set directories here
thisDir = myDirectoryWithTheData;

% load list of experiments for analysis
stimName = 'VisualALL';

load(fullfile(myDirectoryWithTheData,...
    strcat(stimName,'_','expList_forPowerAnalysis.mat')));

AL = expList_forPowerAnalysis; clear expList_forPowerAnalysis

unit = 'na'; %'dB'
minBL = 0.7;        % minimum baseline length in seconds

%% loop through list
for iie = 1:length(AL)
    
    Exps.animal     = AL{iie}.mousename;
    Exps.iseries    = AL{iie}.series;
    Exps.iexp       = AL{iie}.exp;
    
    expRef = strcat(Exps.iseries(1:4),'-',Exps.iseries(5:6),'-',Exps.iseries(7:8),...
        '_',Exps.iexp,'_',Exps.animal);
    
    if length(Exps.iexp) < 2
        [b] = generateGenBlock(expRef, Exps);
    else
        [allExps,~,~,~] = concatenateExps(Exps.animal,expRef(1:10),str2num(Exps.iexp),true,true);
        b = allExps.block;
        t = allExps.t;
        if isfield(b,'exps')      % for concatenated blocks
            b.iexp = num2str(b.exps);
        end
    end
    
    % compute it from visual cortex ROI
    load(fullfile(myDirectoryWithTheData,'ROIPixelSelections',...
        strcat(b.animal,'_',b.iseries,'_',b.iexp,'_PixelPerROI.mat')));
    
    if isfinite(pix(1,1))   % only load experiment if there is an ROI in visual cortex     
        
            %% if power spectra per trial not computed & saved yet
            
            if length(Exps.iexp) < 2        % only one experiment to analyse
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
            
            % changed from 0.7s to 1s on 2017/10/10
            shortbl = find(diff(baseline.OnOffTimes')<minBL);
            
            Fs = round(median(1/median(diff(t))));
            [baseline] = find_noMvmtBaseline(baseline, ntr, Fs);
            
            %% compute power spectra per trial
            
            % make sure chronux toolbox is installed and path is added
            
            pixU = squeeze(U(pix(1,1),pix(2,1),:))';
            avg_trace = pixU*V;
            
            params.tapers = [3 5];
            params.Fs     = round(1/median(diff(t)));
            winSize       = minBL; % in seconds
            
            disp('computing power spectra per trial');
            % compute baseline power spectrum per trial
            for iit = 1:ntr
                if ~logical(sum(ismember(shortbl,iit)))
                    thisTrace = avg_trace(baseline.noMvmtFrames(iit,1):baseline.noMvmtFrames(iit,2));
                    thisTrace = thisTrace - nanmean(thisTrace);
                    if length(thisTrace)*1/params.Fs >= winSize
                        [p, f] = mtspectrumsegc(thisTrace, winSize, params);
                        allPowerSpectra(iit,:) = p';
                    end
                end
            end
            
            allPowerSpectra(allPowerSpectra==0) = NaN;
            if size(allPowerSpectra,1) < ntr
                allPowerSpectra(ntr,:) = NaN;
            end
            
            if ~exist(fullfile(myDirectoryWithTheData,'BaselinePowerSpectra'),'dir');
                mkdir(fullfile(myDirectoryWithTheData,'BaselinePowerSpectra'));
            end
            save(fullfile(myDirectoryWithTheData,'BaselinePowerSpectra',...
                strcat(b.animal,'_',b.iseries,'_',b.iexp,'_Powerspectra_noMvmt_perTrial.mat')),'allPowerSpectra','f');
            
        %%
        
        [~,maxInd] = min(abs(f-10));
        
        [rGO, rNOGO, rGoCorr, rGoInc, rCNG, rGOc, rInZC] = getChoiceNeglectInds(b,ntr,stimName);
        
        PSgo    = allPowerSpectra(rGOc,:);
        PSnogo  = allPowerSpectra(rNOGO,:);
        
        switch unit
            case 'dB'
                
                meanPSgo    = log10(nanmean(PSgo,1));
                meanPSnogo  = log10(nanmean(PSnogo,1));
                
                PSratio     = meanPSgo(1:maxInd) - meanPSnogo(1:maxInd);
                fss = linspace(f(1),f(maxInd),19);
                PSratio_ss = interp1(f(1:maxInd),PSratio,fss);
                
                allRatios(iie,:) = 10*PSratio_ss;   % multiply by 10 to get dB
                
            case 'na'
                
                meanPSgo    = nanmean(PSgo,1);
                meanPSnogo  = nanmean(PSnogo,1);
                
                PSratio     = meanPSgo(1:maxInd) ./ meanPSnogo(1:maxInd);
                fss = linspace(f(1),f(maxInd),19);
                PSratio_ss = interp1(f(1:maxInd),PSratio,fss);
                
                allRatios(iie,:) = PSratio_ss;
                
        end
        
        clear allPowerSpectra allPS PSgo PSnogo meanPSgo meanPSnogo PSratio PSratio_ss b t
        
    end
    switch unit
        case 'dB'
            save(fullfile(myDirectoryWithTheData,...
                strcat(stimName,'_PowerRatio_noMvmt_VisualCortex.mat')),'allRatios','fss');
        case 'na'
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
xticks([0 3 6 9]);
ylim([0.82 1]); 
yticks([0.85 0.9 0.95 1]);
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
