% script with code that generated the behavioural example panels in Figure 1

% visualBehaviour_example

% set directories here
thisDir = myDirectoryWithTheData;

% example dataset used
Exps.animal     = 'Cori';
Exps.iseries    = '20161208';
Exps.iexp       = '1';

expRef = strcat(Exps.iseries(1:4),'-',Exps.iseries(5:6),'-',Exps.iseries(7:8),...
    '_',Exps.iexp,'_',Exps.animal);

[b] = generateGenBlock(expRef, Exps);
[~, ~, t, ~]  = loadSVDfiles(b, true, true);

[baseline, ~] = get_baseline(b, t, false);
shortbl = find(diff(baseline.OnOffTimes')<1);

ntr = b.completedTrials;
if b.excludeFirstTrial
    ntr = ntr-1;
end
tEt   = [b.evts.endTrialTimes(1:ntr)]./60;
tSt   = [b.evts.newTrialTimes(1:ntr)]./60;

Fs  = round(1/median(diff(t)));

frb = [3 6];        % frequency band of interest

nSV = size(U,3);

%% psychometric curves

exampleVisualBehaviour_psychCurve;


%% get behavioural performance across session

bw = 10;        % set sliding window length
[percCorrect, percChoice, percIncorrect, percNogo, percIncorrNogo] = choicesOverTime(b,bw);

xlv = (b.evts.endTrialTimes(end)+1)./60;

figure;
plot([b.evts.newTrialTimes(147)/60 (b.evts.newTrialTimes(147))/60],[0 100],'k');
plot([(b.evts.endTrialTimes(147)+0.5)/60 (b.evts.endTrialTimes(147)+0.5)/60],[0 100],'k');
plot([b.evts.newTrialTimes(285)/60 (b.evts.newTrialTimes(285))/60],[0 100],'k');
plot([(b.evts.endTrialTimes(285)+0.5)/60 (b.evts.endTrialTimes(285)+0.5)/60],[0 100],'k');
bar(tEt,percIncorrNogo,'FaceColor',[0.3 0.3 0.3]); %'k');
plot(tEt,percIncorrect,'Color',[0.64 0.08 0.18],'LineWidth',2);
plot(tEt,percCorrect,'Color',[0 0.9 0.4],'LineWidth',4);
ylabel('Percent');
xlabel('Time (min)');
xlim([0 xlv]);
