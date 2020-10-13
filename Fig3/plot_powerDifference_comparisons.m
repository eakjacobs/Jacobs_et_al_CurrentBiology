% script that plots the Correct - Miss and Incorrect - Miss power
% difference, and Choice - Miss and Correct - Incorrect power difference
% comparisons

% written by Elina Jacobs, UCL Cortexlab

%% initialisations

% set directories here
thisDir = myDirectoryWithTheData;

whichMice   = 'all';
stimName    = 'VisualALL';

% load power differences
for c = 1:4
    switch c
        case 1
            whichCond = 'Choice_Miss';
        case 2
            whichCond = 'Correct_Incorrect';
        case 3
            whichCond = 'Correct_Miss';
        case 4
            whichCond = 'Incorrect_Miss';
    end
    load(fullfile(myDirectoryWithTheData,...
        strcat(stimName,'_',num2str(frb(1)),'to',num2str(frb(2)),'Hz_Power_perROI_',whichCond,'.mat')));
    pDiff = 10*squeeze(-diff(PPP,1,1))';
    switch c
        case 1
            pDiff_ChoiceMiss        = pDiff;
        case 2
            pDiff_CorrectIncorrect  = pDiff;
        case 3
            pDiff_CorrectMiss       = pDiff;
        case 4
            pDiff_IncorrectMiss     = pDiff;
    end
end

%% figure

figure;
suptitle([whichMice ' mice']);
subplot(1,2,1);
plot([-6 5],[-6 5 ],'k--'); hold on;
plot([0 0],[-6 6],':','Color',[0.7 0.7 0.7]);
plot([-6 6],[0 0],':','Color',[0.7 0.7 0.7]);
subplot(1,2,2);
plot([-6 5],[-6 5 ],'k--'); hold on;
plot([0 0],[-6 6],':','Color',[0.7 0.7 0.7]);
plot([-6 6],[0 0],':','Color',[0.7 0.7 0.7]);

ChoiceMiss_avg = nanmean(pDiff_ChoiceMiss,2);
ChoiceMiss_sem = nansem(pDiff_ChoiceMiss,[],2);

CorrectIncorrect_avg = nanmean(pDiff_CorrectIncorrect,2);
CorrectIncorrect_sem = nansem(pDiff_CorrectIncorrect,[],2);

% plot Choice-Miss vs Correct-Incorrect
subplot(1,2,1)
for rr = 1:size(PPP,2)
    switch rr
        case 1
            sym = 'o';
        case 2
            sym = 's';
        case 3
            sym = '<';
        case 4
            sym = 'd';
        case 5
            sym = '>';
    end
    errorbarxy(CorrectIncorrect_avg(rr),ChoiceMiss_avg(rr),...
        CorrectIncorrect_sem(rr),ChoiceMiss_sem(rr),...
        'Color',col,'Marker',sym,'MarkerFaceColor','w','LineWidth',1,'MarkerSize',10);
end
xlabel('Correct - Incorrect');
ylabel('Choice - Miss');
axis square; box off
xlim([-5 1.5]); ylim([-6.1 0.5]);
set(gca,'FontSize',24);

CorrectMiss_avg = nanmean(pDiff_CorrectMiss,2);
CorrectMiss_sem = nansem(pDiff_CorrectMiss,[],2);

IncorrectMiss_avg = nanmean(pDiff_IncorrectMiss,2);
IncorrectMiss_sem = nansem(pDiff_IncorrectMiss,[],2);

% plot Incorrect-Miss vs Correct-Miss
subplot(1,2,2)
for rr = 1:size(PPP,2)
    switch rr
        case 1
            sym = 'o';
        case 2
            sym = 's';
        case 3
            sym = '<';
        case 4
            sym = 'd';
        case 5
            sym = '>';
    end
    errorbarxy(IncorrectMiss_avg(rr),CorrectMiss_avg(rr),...
        IncorrectMiss_sem(rr),CorrectMiss_sem(rr),...
        'Color',col,'Marker',sym,'MarkerFaceColor','w','LineWidth',1,'MarkerSize',10);
end
xlabel('Incorrect - Miss');
ylabel('Correct - Miss');
axis square; box off
xlim([-5 1.5]); ylim([-6.1 0.5]);
set(gca,'FontSize',24);





