% script that creates table in which for each power difference value, there
% is an associated subject, session, genotype & area name so can run a
% statistical model on these factors

% written by Elina Jacobs and Kenneth D. Harris, UCL Cortexlab

%% initialisations

% set directories here
thisDir = myDirectoryWithTheData;

stimName = 'VisualALL';

whichCond   = 'Choice_Miss';
% 'Choice_Miss' 'CorrectReject_Miss' 'Correct_Incorrect' 'FalseAlarm_CorrectReject'
% 'Correct_Miss' 'Incorrect_Miss'
%
% 'Correct_Incorrect_Reward' 'CorrectReject_Incorrect_Reward'

% 'ChoiceMiss_intercDiff' 'ChoiceCorrectReject_intercDiff' 'CorrectRejectMiss_intercDiff'
% 'CorrIncorr_Reward' 'CorrectRejectIncorr_Reward'

doANOVA1 = false;    % whether or not to do one-way ANOVA on pDiff

whichMice = 'all';     % all EJ
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
        switch stimName
            case 'VisualALL'
                expInds = [1:37];
            case 'Visual2ANFC'
                expInds = [1:10];
        end
    case 'NS'
        switch stimName
            case 'VisualALL'
                expInds = [38:58];
            case 'Visual2ANFC'
                expInds = [11:31];
        end
    case 'GCamp6f'
        switch stimName
            case 'VisualALL'
                expInds = [2:19,23:40,57:58];
        end
    case 'GCamp6s'
        switch stimName
            case 'VisualALL'
                expInds = [1,20:22,41:56];
        end
    case '2AFC'
        switch stimName
            case 'VisualALL'
                expInds = [20:37];
        end
    case 'goodPerf'
        switch stimName
            case 'VisualALL'
                expInds = [1:15,17:19,23:25,29:30,33,35,38:44,46:58];
        end
end

% load experiment list
load(fullfile(EJDirs.Data.SVDAResults, 'SUMMARIES','June2017',...
    strcat(stimName,'_','expList_forPowerAnalysis.mat')));
AL = expList_forPowerAnalysis; clear expList_forPowerAnalysis

% load power differences
load(fullfile(myDirectoryWithTheData,...
    strcat(stimName,'_',num2str(frb(1)),'to',num2str(frb(2)),'Hz_Power_perROI_',whichCond,'.mat')));

roiNames = {'VIS','AUD','SS','RSP','MO'};
nROIs = size(PPP,2);

%%

if doANOVA1
    
    pDiff = 10*squeeze(-diff(PPP,1,1))';
    
    cc = 0;         % set counter to zero
    ss = 0;
    for iie = 1:size(pDiff,1)
        thisMouse = AL{iie}.mousename;
        switch thisMouse
            case 'EJ006'
                indicator = 's';
                genotype  = 'Rasgrf_3tg_6s';
            case 'EJ007'
                indicator  = 'f';
                genotype   = 'EMX_3tg_6f';
            case 'EJ009'
                indicator  = 'f';
                genotype   = 'EMX_3tg_6f';
            case 'EJ011';
                indicator  = 'f';
                genotype   = 'EMX_3tg_6f';
            case 'EJ012'
                indicator  = 'f';
                genotype   = 'EMX_3tg_6f';
            case 'EJ013'
                indicator  = 'f';
                genotype   = 'EMX_3tg_6f';
            case 'EJ015'
                indicator  = 'f';
                genotype   = 'EMX_3tg_6f';
            case 'FR053'
                indicator  = 'f';
                genotype   = 'EMX_3tg_6f';
            case 'Muller'
                % Vglut1 (G6f)
                indicator  = 'f';
                genotype   = 'Vglut1_6f';
            case 'Theiler'
                % Vglut1 (G6f)
                indicator  = 'f';
                genotype   = 'Vglut1_6f';
            case 'EJ010'
                % 3tgs GCamp6s
                indicator  = 's';
                genotype   = 'EMX_3tg_6s';
            case 'Cori'
                % tetO G6s
                indicator  = 's';
                genotype   = 'tetO_6s';
            case 'Hench'
                % tetO G6s
                indicator  = 's';
                genotype   = 'tetO_6s';
            case 'Reichstein'
                % tetO G6s
                indicator  = 's';
                genotype   = 'tetO_6s';
            case 'Chain'
                % Snap25 G6s
                indicator  = 's';
                genotype   = 'Snap25_6s';
            case 'Radnitz'
                % Snap25 G6s
                indicator  = 's';
                genotype   = 'Snap25_6s';
        end
        ss = ss+1;
        for iir = 1:nROIs
            if isfinite(pDiff(iie,iir))
                cc = cc+1;
                powerTable(cc).Subject = thisMouse;
                powerTable(cc).Session = ss;
                powerTable(cc).area    = roiNames{iir};
                powerTable(cc).power   = pDiff(iie,iir);
                powerTable(cc).indicator = indicator;
                powerTable(cc).genotype  = genotype;
            end
        end
    end
    
    
    % test ANOVA on power difference
    xx = [powerTable.power];
    b1 = {powerTable.Subject};
    b2 = [powerTable.Session];
    b3 = {powerTable.area};
    b4 = {powerTable.indicator};
    b5 = {powerTable.genotype};
    
    % mixed effects ANOVA
    nsts = [0 0 0 0;
        0 0 1 0
        0 0 0 1
        0 0 0 0];
    [pval,table,stats] = anovan(xx, {b3,b2,b1,b5}, ...
        'model', 'interaction', 'varnames', {'area', 'session', 'mouse', 'genotype'}, ...
        'nested',nsts, ...
        'sstype', 2);
    
    
    [p,anovatab,stats] = anova1(pDiff);
    if p < 0.05
        [mc] = multcompare(stats);
    end
    ANOVA1 = [];
    ANOVA1.pval = p;
    ANOVA1.table = anovatab;
    ANOVA1.stats = stats;
    if p < 0.05
        ANOVA1.multipleComparisons = mc;
    end
    
    clear powerTable
end

%%

conds = strsplit(whichCond,'_');

cc = 0;         % set counter to zero
ss = 0;
for iie = expInds
    thisMouse = AL{iie}.mousename;
    if isfinite(nanmean(ps(1,:,iie))) && isfinite(nanmean(ps(2,:,iie)))
        ss = ss+1;
        
        switch thisMouse
            case 'EJ006'
                indicator = 's';
                genotype  = 'Rasgrf_3tg_6s';
            case 'EJ007'
                indicator  = 'f';
                genotype   = 'EMX_3tg_6f';
            case 'EJ009'
                indicator  = 'f';
                genotype   = 'EMX_3tg_6f';
            case 'EJ011';
                indicator  = 'f';
                genotype   = 'EMX_3tg_6f';
            case 'EJ012'
                indicator  = 'f';
                genotype   = 'EMX_3tg_6f';
            case 'EJ013'
                indicator  = 'f';
                genotype   = 'EMX_3tg_6f';
            case 'EJ015'
                indicator  = 'f';
                genotype   = 'EMX_3tg_6f';
            case 'FR053'
                indicator  = 'f';
                genotype   = 'EMX_3tg_6f';
            case 'Muller'
                % Vglut1 (G6f)
                indicator  = 'f';
                genotype   = 'Vglut1_6f';
            case 'Theiler'
                % Vglut1 (G6f)
                indicator  = 'f';
                genotype   = 'Vglut1_6f';
            case 'EJ010'
                % 3tgs GCamp6s
                indicator  = 's';
                genotype   = 'EMX_3tg_6s';
            case 'Cori'
                % tetO G6s
                indicator  = 's';
                genotype   = 'tetO_6s';
            case 'Hench'
                % tetO G6s
                indicator  = 's';
                genotype   = 'tetO_6s';
            case 'Reichstein'
                % tetO G6s
                indicator  = 's';
                genotype   = 'tetO_6s';
            case 'Chain'
                % Snap25 G6s
                indicator  = 's';
                genotype   = 'Snap25_6s';
            case 'Radnitz'
                % Snap25 G6s
                indicator  = 's';
                genotype   = 'Snap25_6s';
        end
        for iir = 1:nROIs
            for iic = 1:2
                if isfinite(PPP(iic,iir,iie))
                    cc = cc+1;
                    powerTable(cc).Subject  = thisMouse;
                    powerTable(cc).Session  = ss;
                    powerTable(cc).area     = roiNames{iir};
                    powerTable(cc).cond     = conds{iic};
                    powerTable(cc).power    = PPP(iic,iir,iie);
                    powerTable(cc).indicator    = indicator;
                    powerTable(cc).genotype     = genotype;
                end
            end
        end
    end
end


tab = table([powerTable.power]', {powerTable.Subject}', [powerTable.Session]', ...
    {powerTable.area}', {powerTable.cond}', {powerTable.genotype}', ...
    'VariableNames', {'power', 'subject', 'session', 'area', 'cond', 'genotype'});

lme = fitlme(tab, ...
    'power ~ cond + area +  (area|genotype) + (area|subject:genotype) + (area|session:subject:genotype) ');
%, 'CovariancePattern', {'Isotropic','Isotropic','Isotropic'})

anova(lme)

lme = fitlme(tab, ...
    'power ~ cond * area +  (area|genotype) + (area|subject:genotype) + (area|session:subject:genotype) ');
%, 'CovariancePattern', {'Isotropic','Isotropic','Isotropic'})

anova(lme)

