% MissSequenceAnalysis_forSummary

% set directories here
thisDir = myDirectoryWithTheData;

% load experiment list
stimName = 'VisualALL';
load(fullfile(myDirectoryWithTheData,'experimentLists',...
    strcat(stimName,'_','expList_forPowerAnalysis.mat')));
AL = expList_forPowerAnalysis; clear expList_forPowerAnalysis

% loop through experiments
for iie = 1:length(AL) % 1:length(AL) % [1:5,7:length(AL)] for visual expList
    
    Exps.animal     = AL{iie}.mousename;
    Exps.iseries    = AL{iie}.series;
    Exps.iexp       = AL{iie}.exp;
    
    expRef = strcat(Exps.iseries(1:4),'-',Exps.iseries(5:6),'-',Exps.iseries(7:8),...
        '_',Exps.iexp,'_',Exps.animal);
    
    if length(Exps.iexp) < 2        % only one experiment to analyse
        [b] = generateGenBlock(expRef, Exps);
        ntr = b.completedTrials;
        if b.excludeFirstTrial
            ntr = ntr-1;
        end
        
    else
        [allExps,~,~,~] = concatenateExps(Exps.animal,expRef(1:10),str2num(Exps.iexp),true,true);
        b = allExps.block;
        ntr = b.completedTrials;
    end
    
    ee = [b.evts];
    
    if isfield(b,'excludeFirstTrial');
        if b.excludeFirstTrial
            rV = [ee.responseValues(2:ntr)];
            fV = [ee.feedbackValues(2:ntr)];
        else
            rV = [ee.responseValues(1:ntr)];
            fV = [ee.feedbackValues(1:ntr)];
        end
    else
        rV = [ee.responseValues(1:ntr)];
        fV = [ee.feedbackValues(1:ntr)];
    end
    
    rGO  = find(rV~=0);
    rNoGo = find(rV==0);
    fNoGo = fV(rNoGo);          % gives what feedback animal received for giving NoGo response
    fNoGoPos = find(fNoGo>0);
    fNoGoNeg = find(fNoGo<0);
    rNOGO = rNoGo(fNoGoNeg);
    

    %% compute Miss blocks
    
    thr = 1;       % threshold of how many consecutive trials ought to be nogo to be counted as a block
    
    xx = diff(rNOGO);
    gaps = find(xx>thr);          
    % this gives vector with indices with rNOGO that are the last trials of
    % a NOGO block
    singleNOGOs = find(xx==1);
    
    nogoBlockLengths = [];
    cc = 0;
    if rNOGO(1) == 1
        if rNOGO(gaps(1)) > 1
            cc = cc+1;
            nogoBlockLengths(cc) = rNOGO(gaps(1));
        end
    else 
        thisBlockStart = rNOGO(1);
        thisBlockEnd   = rNOGO(gaps(1));
        if thisBlockEnd - thisBlockStart > 0
            if thisBlockEnd - thisBlockStart + 1 > 1
                cc = cc+1;
                nogoBlockLengths(cc) = thisBlockEnd - thisBlockStart + 1; % +1 since want to include beginning and end
            end
        end
    end
    for iin = 2:length(gaps)
        thisBlockStart = rNOGO(gaps(iin-1)+1);
        thisBlockEnd   = rNOGO(gaps(iin));
        if thisBlockEnd - thisBlockStart > 0
            cc = cc+1;
            nogoBlockLengths(cc) = thisBlockEnd - thisBlockStart + 1; % +1 since want to include beginning and end
        end
    end
    
    nMissPeriods = length(nogoBlockLengths);
    
    %% shuffle trials, compute shuffled Miss blocks
    
    shuffledTrials = randperm(ntr);
    rVsh = rV(shuffledTrials);
    fVsh = fV(shuffledTrials);
    
    rGOsh  = find(rVsh~=0);
    rNoGosh = find(rVsh==0);
    fNoGosh = fVsh(rNoGosh);          % gives what feedback animal received for giving NoGo response
    fNoGoPossh = find(fNoGosh>0);
    fNoGoNegsh = find(fNoGosh<0);
    rNOGOsh = rNoGo(fNoGoNegsh);
    
    
    xSh = diff(rNOGOsh);
    gapsSh = find(xSh>thr);
    % this gives vector with indices with rNOGOsh that are the last trials of
    % a NOGO block
    singleNOGOs_sh = find(xSh==1);
    
    nogoBlockLengths_sh = [];
    cc = 0;
    if rNOGOsh(1) == 1
        if rNOGOsh(gapsSh(1)) > 1
            cc = cc+1;
            nogoBlockLengths_sh(cc) = rNOGOsh(gapsSh(1));
        end
    else 
        thisBlockStart_sh = rNOGOsh(1);
        thisBlockEnd_sh   = rNOGOsh(gapsSh(1));
        if thisBlockEnd_sh - thisBlockStart_sh > 0
            if thisBlockEnd_sh - thisBlockStart_sh + 1 > 1
                cc = cc+1;
                nogoBlockLengths_sh(cc) = thisBlockEnd_sh - thisBlockStart_sh + 1; % +1 since want to include beginning and end
            end
        end
    end
    for iin = 2:length(gapsSh)
        thisBlockStart_sh = rNOGOsh(gapsSh(iin-1)+1);
        thisBlockEnd_sh   = rNOGOsh(gapsSh(iin));
        if thisBlockEnd_sh - thisBlockStart_sh > 0
            cc = cc+1;
            nogoBlockLengths_sh(cc) = thisBlockEnd_sh - thisBlockStart_sh + 1; % +1 since want to include beginning and end
        end
    end
    
    nMissPeriods_sh = length(nogoBlockLengths_sh);
    
    
    %% summary
    
    MissSequences(iie).nMissTrials    = length(rNOGO);
    MissSequences(iie).totalTrials       = ntr;
    
    MissSequences(iie).nMissPeriods   = nMissPeriods;
    MissSequences(iie).MissPeriodLengths = nogoBlockLengths;
    MissSequences(iie).nSingleNogos      = length(singleNOGOs);
    
    MissSequences(iie).Shuffled_nMissPeriods   = nMissPeriods_sh;
    MissSequences(iie).Shuffled_MissPeriodLengths = nogoBlockLengths_sh;
    MissSequences(iie).Shuffled_nSingleNogos      = length(singleNOGOs_sh);
    
    clear xx xSh gaps gapsSh rNOGO rNOGOsh nogoBlockLengths nogoBlockLengths_sh nMissPeriods nMissPeriods_sh thisBlockStart thisBlockStart_sh thisBlockEnd thisBlockEnd_sh
    
    
end

%% figures

% scatterplot
figure; 
plot([0 450],[0 450],':','Color',[0.7 0.7 0.7]);  hold on
text(301,308,'y=x','FontSize',14);
plot([0 450],[0 450/2],':','Color',[0.7 0.7 0.7]);
text(453, 450/2,'y=x/2','FontSize',14);
plot([0 450],[0 450/4],':','Color',[0.7 0.7 0.7]);
text(453, 450/4,'y=x/4','FontSize',14);
plot([neglectSequences.totalTrials],[neglectSequences.nNeglectTrials],'ko','MarkerSize',48);
axis square; box off
xlim([0 450]); ylim([0 300]);
xlabel('Total trials');
ylabel('Neglect trials');
set(gca, 'FontSize', 24)
x0=200;
y0=200;
width = 300;
height= 300;
set(gca,'units','points','position',[x0,y0,width,height]);

% histograms
figure;
subplot(2,3,1)
h1 = histogram([MissSequences.nMissPeriods],30);
title('Number of Miss Periods per dataset');
box off

subplot(2,3,2)
h2 = histogram([MissSequences.MissPeriodLengths],129);
box off
title('Miss Period Lengths');

subplot(2,3,3)
h3 = histogram([MissSequences.nSingleNogos],50);
box off
title('Single Nogos');

subplot(2,3,4)
h4 = histogram([MissSequences.Shuffled_nMissPeriods],30);
title('Number of Miss Periods after shuffling trials');
box off

subplot(2,3,5)
h5 = histogram([MissSequences.Shuffled_MissPeriodLengths],129);
box off
title('Shuffled Miss Period Lengths');

subplot(2,3,6)
h3 = histogram([MissSequences.Shuffled_nSingleNogos],50);
box off
title('Shuffled Single Nogos');