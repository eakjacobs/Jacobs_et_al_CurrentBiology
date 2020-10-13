% exampleVisualBehaviour_psychCurve

% plots the psychometric curve for the 2AFC trials of the example in
% behavioural figure

% written by Elina Jacobs, UCL Cortexlab
% this script requires the NeuroGLM repository, 
% written by Peter Zatka-Haas: https://github.com/peterzh/NeuroGLM

[EJDirs] = setEJDirs;
% SetDefaultDirs;
% BehDir   = DIRS.expInfo;
BehDir = '\\zclone.cortexlab.net\Data\expInfo';

Exps.animal     = 'Cori';
Exps.iseries    = '20161208';
Exps.iexp       = '1';


expRef = strcat(Exps.iseries(1:4),'-',Exps.iseries(5:6),'-',Exps.iseries(7:8),...
    '_',Exps.iexp,'_',Exps.animal);

[b] = generateGenBlock(expRef, Exps);

%% get variables for model
% get rid of trials with contrasts on both sides
zeroCleft   = find(b.stimuli(:,1)==0);
zeroCright  = find(b.stimuli(:,2)==0);
AFC         = unique(sort([zeroCleft;zeroCright]));     % gives indices of trials in which at least one of contrasts was zero
while AFC(end)>b.completedTrials
    AFC(end) = [];
end
stimuli     = b.stimuli(AFC,:);

response            = [b.evts.responseValues(AFC)];        % -1 when turned right, ie acted like there was a leftvis/highfreq stim, 1 when turned left ie acted like there was rigtvis/low freq stim
repeatNum           = [b.evts.repeatNum(AFC)];       % 1 when first trial, 2 when second, and so on
repeatTrials        = find(repeatNum>1);

conditions          = diff(stimuli');
stimLeft            = find(conditions<0);
stimRight           = find(conditions>0);

firstAttempts               = response; % -1 when turned right (ie like high stim), 1 when turned left (ie like low stim)
firstAttempts(repeatTrials) = NaN;
firstAttempts               = firstAttempts(isfinite(firstAttempts));
% now put it into format that it can go into D for Peter ZH's model
r = firstAttempts';
% uncomment to use Peter ZH's model
r(firstAttempts==0) = 3;
r(firstAttempts==1) = 2;
r(firstAttempts==-1)= 1;

cond                = conditions;
cond(repeatTrials)  = NaN;
cond                = cond(isfinite(cond));
% now put it into format that it can go into D for Peter ZH's model
c = zeros(length(cond),2);
cl = find(cond<0); cr = find(cond>0);
c(cl,1) = abs(cond(cl));
c(cr,2) = cond(cr);

f = [b.evts.feedbackValues(AFC)]';
f(repeatTrials) = NaN;
f = f(isfinite(f));

rN = [b.evts.repeatNum(AFC)]';
rN(repeatTrials) = NaN; rN = rN(isfinite(rN));
rN(rN==1) = 0;
rN(rN>1) = 1;

D.contrast_cond = c;
D.response      = r;
D.feedbackType  = f;
D.repeatNum     = rN;
% to get errorbars using EJ's function
stim_cond_EJ    = cond;
response_EJ     = firstAttempts';


%% run model & plot

unique_conditions  = unique(stim_cond_EJ);
cLeft = find(unique_conditions<0);
cZero = find(unique_conditions==0);
cRight= find(unique_conditions>0);
ucl = length(unique_conditions);

for xx = 1:ucl
    if unique_conditions(xx) < 0
        stimSide = 'L';
    elseif unique_conditions(xx) > 0
        stimSide = 'R';
    else
        stimSide = '';
    end
    xLabz{xx} = strcat(stimSide,num2str(abs(unique_conditions(xx))*100));
end

[~,~,~,percent_right,percent_left,percent_NoGo,ste_right,ste_left,ste_NoGo,~,~,~] = ...
    psych_curve_nogo(response_EJ,stim_cond_EJ,unique_conditions,ucl);

figure;
plot([0 0],[0 1],'k:'); hold on;
plot([-1 1],[0.5 0.5],'k:'); hold on;
errorbar(unique_conditions, percent_left, ste_left, 'ko', 'markerfacec','k'); hold on;
errorbar(unique_conditions, percent_right, ste_right, 'ko', 'markerfacec','k'); hold on;
errorbar(unique_conditions, percent_NoGo, ste_right, 'ko', 'markerfacec','k'); hold on;

gv = GLM(D);
gv = gv.setModel('C50-subset');
gv = gv.fit;
gv.plotFit; hold on;

xlim([-1 1]);
ylim([0 1]);
xticks([-1 -0.5 -0.25 0 0.25 0.5 1])
xticklabels({'L100','L50','L25','0','R25','R50','R100'});

set(gca, 'FontSize', 18);
axis square;

%%

maxC = max(max(gv.data.stimulus));
evalC = [linspace(maxC,0,100)', zeros(100,1);
    zeros(100,1), linspace(0,maxC,100)'];
evalC1d = evalC(:,2) - evalC(:,1);

phat = gv.calculatePhat(gv.parameterFits,evalC);
phatOrdered = zeros(size(phat));
phatOrdered(:,1) = phat(:,2);
phatOrdered(:,2) = phat(:,3);
phatOrdered(:,3) = phat(:,1);

spOrder = [3,2,1];
figure;
suptitle([b.animal ' ' b.iseries ' ' b.iexp]);

for sp = 1:3
    subplot(1,3,spOrder(sp))
    
    plot([0 0],[0 1],'k:'); hold on;
    plot([-1 1],[0.5 0.5],'k:');
    
    switch sp
        case 1
            errorbar(unique_conditions(cRight), percent_left(cRight), ste_left(cRight),...
                'ko', 'markerfacec',[0 0.9 0.4],'MarkerSize',15);
            errorbar(unique_conditions(cZero), percent_left(cZero), ste_left(cZero),...
                'ko', 'markerfacec',[0.64 0.08 0.18],'MarkerSize',15);
            errorbar(unique_conditions(cLeft), percent_left(cLeft), ste_left(cLeft),...
                'ko', 'markerfacec',[0.64 0.08 0.18],'MarkerSize',15);
        case 3
            errorbar(unique_conditions(cRight), percent_right(cRight), ste_right(cRight),...
                'ko', 'markerfacec',[0.64 0.08 0.18],'MarkerSize',15);
            errorbar(unique_conditions(cZero), percent_right(cZero), ste_right(cZero),...
                'ko', 'markerfacec',[0.64 0.08 0.18],'MarkerSize',15);
            errorbar(unique_conditions(cLeft), percent_right(cLeft), ste_right(cLeft),...
                'ko', 'markerfacec',[0 0.9 0.4],'MarkerSize',15);
        case 2
            errorbar(unique_conditions(cRight), percent_NoGo(cRight), ste_right(cRight),...
                'ko', 'markerfacec',[0.3 0.3 0.3],'MarkerSize',15);
            errorbar(unique_conditions(cZero), percent_NoGo(cZero), ste_right(cZero),...
                'ko', 'markerfacec',[0 0.9 0.4],'MarkerSize',15);
            errorbar(unique_conditions(cLeft), percent_NoGo(cLeft), ste_right(cLeft),...
                'ko', 'markerfacec',[0.3 0.3 0.3],'MarkerSize',15);
    end
    
    %     gv.plotFit; hold on;
    plot(evalC1d,phatOrdered(:,sp),'Color',[0.5 0.5 0.5],'LineWidth',1);
    
    xlim([-1.05 1.05]);
    ylim([0 1]);
    xticks([unique_conditions])
    xticklabels(xLabz);
    xlabel('Contrast (%)');
    switch sp
        case 1
            ylabel('P(Right Choice)');
        case 2
            ylabel('P(NoGo)');
        case 3
            ylabel('P(Left Choice)');
    end
    set(gca, 'FontSize', 26);
    axis square; box off;
end

