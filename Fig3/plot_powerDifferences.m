% script that plots the results from compute_powerDifferences

% written by Elina Jacobs, UCL Cortexlab

%% initialisations

% set directories here
thisDir = myDirectoryWithTheData;

pr = 'ROIs';
frb = [3 6];
whichCond = 'Choice_Miss';
stimName = 'VisualALL';
% what experiment type: Visual, Audio or AudioVisual

% load experiment list
load(fullfile(myDirectoryWithTheData,...
    strcat(stimName,'_','expList_forPowerAnalysis.mat')));
AL = expList_forPowerAnalysis; clear expList_forPowerAnalysis

% load results
load(fullfile(myDirectoryWithTheData,...
    strcat(stimName,'_',num2str(frb(1)),'to',num2str(frb(2)),'Hz_Power_perROI_',whichCond,'.mat')));


whichMice = 'all';     % can sub select datasets according to criteria; only applicable for visual datasets
switch whichMice
    case 'all'
        switch stimName
            case 'VisualALL'
                expInds = [1:58];
            case 'Visual2ANFC'
                expInds = [1:31];
            case 'Auditory'
                expInds = [1:15];
            case 'AudioVisual'
                expInds = [1:10];
        end
    case 'EJ'
        expInds = [1:37];
    case 'NS'
        expInds = [38:58];
    case 'GCamp6s'
        expInds = [1,20:22,41:56];
    case 'GCamp6f'
        expInds = [2:19,23:40,57:58];
    case 'goodPerf'         % datasets in which the percent incorrect at maximum contrast was less than 10%
        expInds = [1:15,17:19,23:25,29:30,33,35,38:44,46:58];
    case '2AFC'
        expInds = [20:37];
end


% need to specify ROIs and order for plotting purposes
switch stimName       
        
    case 'Auditory'
        roiNames{2} = 'VIS';
        roiNames{3} = 'AUD';
        roiNames{1} = 'SS';
        roiNames{4} = 'RSP';
        
        ROIorder = [3 1 2 4];
        
        
    case 'AudioVisual'
        roiNames{2} = 'VIS';
        roiNames{3} = 'AUD';
        roiNames{1} = 'SS';
        roiNames{4} = 'RSP';
        
        ROIorder = [3 1 2 4];
        
        
    case 'Visual2ANFC'
        roiNames{2} = 'VIS';
        roiNames{4} = 'AUD';
        roiNames{1} = 'SS';
        roiNames{5} = 'RSP';
        roiNames{3} = 'MO';
        
        ROIorder = [3 1 5 2 4];
        
        
    case 'VisualALL'
        
        switch whichMice
            case 'NS'
                roiNames{2} = 'VIS';
                roiNames{4} = 'AUD';
                roiNames{1} = 'SS';
                roiNames{5} = 'RSP';
                roiNames{3} = 'MO';
                
                ROIorder = [3 1 5 4];
                
            case 'NS_noHench'
                roiNames{2} = 'VIS';
                roiNames{4} = 'AUD';
                roiNames{1} = 'SS';
                roiNames{5} = 'RSP';
                roiNames{3} = 'MO';
                
                ROIorder = [3 1 5 4];
                
            otherwise
                roiNames{2} = 'VIS';
                roiNames{4} = 'AUD';
                roiNames{1} = 'SS';
                roiNames{5} = 'RSP';
                roiNames{3} = 'MO';
                
                ROIorder = [3 1 5 2 4];
        end
        
end

nROIs = length(ROIorder);

%% 

switch pr
    case 'ROIs'
        
        pDiff = 10*squeeze(-diff(PPP,1,1))';
        
        fs = figure;
        for iir = 1:nROIs; % [3 1 5 2 4];  % to plot ROIs in order SS - Vis - MOS - Aud - RS
            rInd = ROIorder(iir);
            allExpsPD = pDiff(:,rInd);
            thisExpPD = allExpsPD(expInds);
            
            % all symbols: +o*.xsd^v<>ph
            for iie = expInds
                switch AL{iie}.mousename
                    case 'EJ006'
                        xp = 0;
                        ms = '\heartsuit';
                        text([xp+2*iir],allExpsPD(iie),ms,'Color','b','FontSize',12);
                    case 'EJ007'
                        ms = '+';
                        plot([0.07+2*iir],allExpsPD(iie),ms,'Color',[0.85 0.85 0.85],'MarkerSize',10);
                    case 'EJ009'
                        ms = 'o';
                        plot([0.14+2*iir],allExpsPD(iie),ms,'Color',[0.85 0.85 0.85],'MarkerSize',10);
                    case 'EJ011';
                        ms = '*';
                        plot([0.21+2*iir],allExpsPD(iie),ms,'Color',[0.85 0.85 0.85],'MarkerSize',10);
                    case 'EJ012'
                        ms = '.';
                        plot([0.28+2*iir],allExpsPD(iie),ms,'MarkerSize',14,'Color',[0.85 0.85 0.85]);
                    case 'EJ013'
                        ms = 'x';
                        plot([0.35+2*iir],allExpsPD(iie),ms,'Color',[0.85 0.85 0.85],'MarkerSize',10);
                    case 'EJ015'
                        ms = 's';
                        plot([0.42+2*iir],allExpsPD(iie),ms,'Color',[0.85 0.85 0.85],'MarkerSize',10);
                    case 'FR053'
                        ms = 'd';
                        plot([0.49+2*iir],allExpsPD(iie),ms,'Color',[0.85 0.85 0.85],'MarkerSize',10);
                    case 'Muller'
                        % Vglut1 (G6f)
                        ms = '^';
                        plot([0.56+2*iir],allExpsPD(iie),ms,'Color',[0.85 0.85 0.85],'MarkerSize',10);
                    case 'Theiler'
                        % Vglut1 (G6f)
                        ms = 'v';
                        plot([0.63+2*iir],allExpsPD(iie),ms,'Color',[0.85 0.85 0.85],'MarkerSize',10);
                    case 'EJ010'
                        % 3tgs GCamp6s
                        ms = '<';
                        plot([0.7+2*iir],allExpsPD(iie),ms,'Color',[0.85 0.85 0.85],'MarkerSize',10);
                    case 'Cori'
                        % tetO G6s
                        ms = '>';
                        plot([0.77+2*iir],allExpsPD(iie),ms,'Color',[0.85 0.85 0.85],'MarkerSize',10);
                    case 'Hench'
                        % tetO G6s
                        ms = 'p';
                        plot([0.84+2*iir],allExpsPD(iie),ms,'Color',[0.85 0.85 0.85],'MarkerSize',10);
                    case 'Reichstein'
                        % tetO G6s
                        ms = 'h';
                        plot([0.91+2*iir],allExpsPD(iie),ms,'Color',[0.85 0.85 0.85],'MarkerSize',10);
                    case 'Chain'
                        % Snap25 G6s
                        xp = 0.98;
                        ms = '\spadesuit';
                        text([xp+2*iir],allExpsPD(iie),ms,'Color',[0.9 0.9 0.9],'FontSize',12);
                    case 'Radnitz'
                        % Snap25 G6s
                        xp = 1.05;
                        ms = '\clubsuit';
                        text([xp+2*iir],allExpsPD(iie),ms,'Color',[0.9 0.9 0.9],'FontSize',12);
                end
                
            end
            
            %             % to test if normally distributed
            %             test_cdf = [thisPD,cdf('tlocationscale',thisPD,nanmean(thisPD),std(thisPD(isfinite(thisPD))),1)];
            %             kstest(thisPD,'CDF',test_cdf)
            
            % plot overall average
            offset = 1.5;
            if iir == 4
                switch stimName
                    case 'Visual2ANFC'
                        offset = 0.6;
                    otherwise
                        offset = 1.1;
                end
                switch whichMice
                    case 'NS'
                        offset = 1.5;
                end
            end
            switch stimName
                case 'Auditory'
                    offset = 0.7;
                case 'AudioVisual'
                    offset = 0.7;
            end
            plot([offset+2*iir],nanmean(thisExpPD),'k.','MarkerSize',18);
            errorbar([offset+2*iir],nanmean(thisExpPD),SEM(thisExpPD(isfinite(thisExpPD))),'Color','k','LineWidth',1);
            
            % plot average per animal
            switch stimName
                case 'VisualALL'
                    
                    switch whichMice
                        case 'all'
                            xp = 0;
                            ms = '\heartsuit';
                            text([xp+2*iir],allExpsPD(1),ms,'Color','b','FontSize',12);
                            
                            EJ7inds     = [2:5,11:15];
                            ms = '+';
                            plot([0.07+2*iir],nanmean(allExpsPD(EJ7inds)),ms,'Color',[0 0.8 0.6],'MarkerSize',10);
                            
                            EJ9inds     = [6:10,16:19];
                            ms = 'o';
                            plot([0.14+2*iir],nanmean(allExpsPD(EJ9inds)),ms,'Color',[0 0.8 0.6],'MarkerSize',10);
                            
                            EJ10inds    = [20:22];
                            ms = '<';
                            plot([0.7+2*iir],nanmean(allExpsPD(EJ10inds)),ms,'Color','c','MarkerSize',10);
                            
                            FRinds      = [23:25];
                            ms = 'd';
                            plot([0.49+2*iir],nanmean(allExpsPD(FRinds)),ms,'Color',[0 0.8 0.6],'MarkerSize',10);
                            
                            EJ12inds    = [26:27];
                            ms = '.';
                            plot([0.28+2*iir],nanmean(allExpsPD(EJ12inds)),ms,'MarkerSize',14,'Color',[0 0.8 0.6]);
                            
                            EJ13inds    = [28:31];
                            ms = 'x';
                            plot([0.35+2*iir],nanmean(allExpsPD(EJ13inds)),ms,'Color',[0 0.8 0.6],'MarkerSize',10);
                            
                            EJ15inds    = [32:37];
                            ms = 's';
                            plot([0.42+2*iir],nanmean(allExpsPD(EJ15inds)),ms,'Color',[0 0.8 0.6],'MarkerSize',10);
                            
                        case '2AFC'
                            
                            EJ10inds    = [20:22];
                            ms = '<';
                            plot([0.7+2*iir],nanmean(allExpsPD(EJ10inds)),ms,'Color','c','MarkerSize',10);
                            
                            FRinds      = [23:25];
                            ms = 'd';
                            plot([0.49+2*iir],nanmean(allExpsPD(FRinds)),ms,'Color',[0 0.8 0.6],'MarkerSize',10);
                            
                            EJ12inds    = [26:27];
                            ms = '.';
                            plot([0.28+2*iir],nanmean(allExpsPD(EJ12inds)),ms,'MarkerSize',14,'Color',[0 0.8 0.6]);
                            
                            EJ13inds    = [28:31];
                            ms = 'x';
                            plot([0.35+2*iir],nanmean(allExpsPD(EJ13inds)),ms,'Color',[0 0.8 0.6],'MarkerSize',10);
                            
                            EJ15inds    = [32:37];
                            ms = 's';
                            plot([0.42+2*iir],nanmean(allExpsPD(EJ15inds)),ms,'Color',[0 0.8 0.6],'MarkerSize',10);
                            
                        case 'goodPerf'
                            xp = 0;
                            ms = '\heartsuit';
                            text([xp+2*iir],allExpsPD(1),ms,'Color','b','FontSize',12);
                            
                            EJ7inds     = [2:5,11:15];
                            ms = '+';
                            plot([0.07+2*iir],nanmean(allExpsPD(EJ7inds)),ms,'Color',[0 0.8 0.6],'MarkerSize',10);
                            
                            EJ9inds     = [6:10,17:19];
                            ms = 'o';
                            plot([0.14+2*iir],nanmean(allExpsPD(EJ9inds)),ms,'Color',[0 0.8 0.6],'MarkerSize',10);
                            
                            FRinds      = [23:25];
                            ms = 'd';
                            plot([0.49+2*iir],nanmean(allExpsPD(FRinds)),ms,'Color',[0 0.8 0.6],'MarkerSize',10);
                            
                            EJ13inds    = [29:30];
                            ms = 'x';
                            plot([0.35+2*iir],nanmean(allExpsPD(EJ13inds)),ms,'Color',[0 0.8 0.6],'MarkerSize',10);
                            
                            EJ15inds    = [33,35];
                            ms = 's';
                            plot([0.42+2*iir],nanmean(allExpsPD(EJ15inds)),ms,'Color',[0 0.8 0.6],'MarkerSize',10);
                    end
                    
                    if ~isequal(whichMice,'2AFC')
                        Mullerinds  = [38:40];
                        ms = '^';
                        plot([0.56+2*iir],nanmean(allExpsPD(Mullerinds)),ms,'Color',[0 0 0.3],'MarkerSize',10);
                        
                        Chaininds   = [41:43];
                        xp = 0.98;
                        ms = '\spadesuit';
                        text([xp+2*iir],nanmean(allExpsPD(Chaininds)),ms,'Color',[0.5 0.18 0.56],'FontSize',12);
                        
                        switch whichMice
                            case 'goodPerf'
                                Coriinds = [44,46:48];
                            otherwise
                                Coriinds = [44:48];
                        end
                        ms = '>';
                        plot([0.77+2*iir],nanmean(allExpsPD(Coriinds)),ms,'Color','m','MarkerSize',10);
                        
                        Radinds     = [49:52];
                        xp = 1.05;
                        ms = '\clubsuit';
                        text([xp+2*iir],nanmean(allExpsPD(Radinds)),ms,'Color',[0.5 0.18 0.56],'FontSize',12);
                        
                        switch whichMice
                            case 'all'
                                Henchinds   = [53:55];
                                ms = 'p';
                                plot([0.84+2*iir],nanmean(allExpsPD(Henchinds)),ms,'Color','m','MarkerSize',10);
                            case 'goodPerf'
                                Henchinds   = [53:55];
                                ms = 'p';
                                plot([0.84+2*iir],nanmean(allExpsPD(Henchinds)),ms,'Color','m','MarkerSize',10);
                            case 'GCamp6s'
                                Henchinds   = [53:55];
                                ms = 'p';
                                plot([0.84+2*iir],nanmean(allExpsPD(Henchinds)),ms,'Color','m','MarkerSize',10);
                                EJ10inds    = [20:22];
                                ms = '<';
                                plot([0.7+2*iir],nanmean(allExpsPD(EJ10inds)),ms,'Color','c','MarkerSize',10);
                            case 'NS'
                                Henchinds   = [53:55];
                                ms = 'p';
                                plot([0.84+2*iir],nanmean(allExpsPD(Henchinds)),ms,'Color','m','MarkerSize',10);
                        end
                        
                        Reichinds   = [56];
                        ms = 'h';
                        plot([0.91+2*iir],allExpsPD(Reichinds),ms,'Color','m','MarkerSize',10);
                        
                        Theilerinds = [57:58];
                        ms = 'v';
                        plot([0.63+2*iir],nanmean(allExpsPD(Theilerinds)),ms,'Color',[0 0 0.3],'MarkerSize',10);
                    end
                    
                case 'Visual2ANFC'
                    
                    xp = 0;
                    ms = '\heartsuit';
                    text([xp+2*iir],allExpsPD(1),ms,'Color','b','FontSize',12);
                    
                    EJ7inds     = [2:5];
                    ms = '+';
                    plot([0.07+2*iir],nanmean(allExpsPD(EJ7inds)),ms,'Color',[0 0.8 0.6],'MarkerSize',10);
                    
                    EJ9inds     = [6:10];
                    ms = 'o';
                    plot([0.14+2*iir],nanmean(allExpsPD(EJ9inds)),ms,'Color',[0 0.8 0.6],'MarkerSize',10);
                    
                    Mullerinds  = [11:13];
                    ms = '^';
                    plot([0.56+2*iir],nanmean(allExpsPD(Mullerinds)),ms,'Color',[0 0 0.3],'MarkerSize',10);
                    
                    Chaininds   = [14:16];
                    xp = 0.98;
                    ms = '\spadesuit';
                    text([xp+2*iir],nanmean(allExpsPD(Chaininds)),ms,'Color',[0.5 0.18 0.56],'FontSize',12);
                    
                    Coriinds    = [17:21];
                    ms = '>';
                    plot([0.77+2*iir],nanmean(allExpsPD(Coriinds)),ms,'Color','m','MarkerSize',10);
                    
                    Reichinds   = [29];
                    ms = 'h';
                    plot([0.91+2*iir],allExpsPD(Reichinds),ms,'Color',[0.85 0.85 0.85],'MarkerSize',10);
                    
                    Radinds     = [22:25];
                    xp = 1.05;
                    ms = '\clubsuit';
                    text([xp+2*iir],nanmean(allExpsPD(Radinds)),ms,'Color',[0.5 0.18 0.56],'FontSize',12);
                    
                    Henchinds   = [26:28];
                    ms = 'p';
                    plot([0.84+2*iir],nanmean(allExpsPD(Henchinds)),ms,'Color','m','MarkerSize',10);
                    
                    Theilerinds = [30:31];
                    ms = 'v';
                    plot([0.63+2*iir],nanmean(allExpsPD(Theilerinds)),ms,'Color',[0 0 0.3],'MarkerSize',10);
                    
                case 'Auditory'
                    
                    EJ11inds    = [1:9];
                    ms = '*';
                    plot([0.21+2*iir],nanmean(allExpsPD(EJ11inds)),ms,'Color',[0 0.8 0.6],'MarkerSize',10);
                    
                    EJ13inds    = [10:14];
                    ms = 'x';
                    plot([0.35+2*iir],nanmean(allExpsPD(EJ13inds)),ms,'Color',[0 0.8 0.6],'MarkerSize',10);
                    
                    EJ15inds    = 15;
                    ms = 's';
                    plot([0.42+2*iir],nanmean(allExpsPD(EJ15inds)),ms,'Color',[0 0.8 0.6],'MarkerSize',10);
                    
                case 'AudioVisual'
                    
                    EJ7inds    = [1:7];
                    ms = '+';
                    plot([0.07+2*iir],nanmean(allExpsPD(EJ7inds)),ms,'Color',[0 0.8 0.6],'MarkerSize',10);
                    
                    EJ9inds    = [8:10];
                    ms = 'o';
                    plot([0.14+2*iir],nanmean(allExpsPD(EJ9inds)),ms,'Color',[0 0.8 0.6],'MarkerSize',10);
                    
                    
            end
        end
        
        xlim([1.5 10]);
        switch stimName
            case 'VisualALL'
                xlim([1.5 12]);
            case 'Visual2ANFC'
                xlim([1.5 12]);
        end
        switch whichMice
            case 'NS'
                xlim([1.5 10]);
        end
        ylim([-8.5 8]);
        switch whichCond
            case 'Choice_Miss'
                switch stimName
                    case 'Auditory'
                        ylim([-10 6.5]);
                    case 'AudioVisual'
                        ylim([-9 7.5]);
                end
            case 'Correct_Incorrect'
                switch stimName
                    case 'AudioVisual'
                        ylim([-11 11]);
                end
        end

        xticks([2.5 4.5 6.5 8.5 10.5]);
        yticks([-5 0 5]);
        box OFF;
        set(gca,'xcolor','none')
        ypos = -9.5;
        text(2.3,ypos,roiNames{1},'FontSize',18);
        text(4.3,ypos,roiNames{2},'FontSize',18);
        text(6.4,ypos,roiNames{3},'FontSize',18);
        switch whichMice
            case 'NS'
                text(8.3,ypos,roiNames{5},'FontSize',18);
            otherwise
                text(8.3,ypos,roiNames{4},'FontSize',18);
                switch stimName
                    case 'VisualALL'
                        text(10.3,ypos,roiNames{5},'FontSize',18);
                    case 'Visual2ANFC'
                        text(10.3,ypos,roiNames{5},'FontSize',18);
                end
        end
        ylabel('Power (dB)');
        conds = strsplit(whichCond,'_');
        title([conds{1} '-' conds{2}]);
        set(gca, 'FontSize', 18)
        x0=400;
        y0=300;
        switch stimName
            case 'Auditory'
                width=350;
            case 'AudioVisual'
                width=350;
            otherwise
                width=400;
        end
        height=180;
        set(gca,'units','points','position',[x0,y0,width,height]);
        
        
        
    case 'average'
        
        pDiff = 10*-diff(PPP);
        [h p] = ttest(pDiff);
        
        fs = figure;
        plot([0 4],[0 0],'k:'); hold on;
        plot([rand(1,length(pDiff))+1],pDiff,ms,'MarkerEdgeColor',[0.5 0.5 0.6]);
        tt = strcat('p=',num2str(p));
        text(1,-9,tt);
        xlim([0 3]);
        ylim([-10.5 5]);
        xticks([1.5]);
        xticklabels('average');
        yticks([-10 -5 0]);
        ylabel('Choice - Neglect 3-6Hz Power');
        title(['Summary across ' stimName ' experiments']);
        
end
