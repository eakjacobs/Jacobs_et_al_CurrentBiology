% script that plots the example auditory psychometric curve

% written by Elina Jacobs, UCL Cortexlab
% this script requires the NeuroGLM repository, 
% written by Peter Zatka-Haas: https://github.com/peterzh/NeuroGLM

% set directories here
thisDir = myDirectoryWithTheData;

load(fullfile(myDirectoryWithTheData,'EJ011_expList_forAnalysis.mat'));

vv = 0;
D = [];
D.contrast_cond     = [];
D.response          = [];
D.repeatNum         = [];
D.feedbackType      = [];

%%

for iie = 4:length(Analysis_ExpList)   % earlier datasets have biased performance
    
    if size(Analysis_ExpList(iie).ExpRef,2) > 18
        Exps.animal  = Analysis_ExpList(iie).ExpRef(15:end);
        Exps.iseries = Analysis_ExpList(iie).ExpRef(1:10);
        Exps.iexp    = Analysis_ExpList(iie).ExpRef(12:13);
    else
        Exps.animal  = Analysis_ExpList(iie).ExpRef(14:end);
        Exps.iseries = Analysis_ExpList(iie).ExpRef(1:10);
        Exps.iexp    = Analysis_ExpList(iie).ExpRef(12);
    end
    
    [b] = generateGenBlock(Analysis_ExpList(iie).ExpRef,Exps);
    
    ntr = b.completedTrials;
    
    response            = [b.evts.responseValues(1:ntr)];       % -1 when turned right, ie acted like there was a leftvis/highfreq stim, 1 when turned left ie acted like there was rigtvis/low freq stim
    repeatNum           = [b.evts.repeatNum(1:ntr)];            % 1 when first trial, 2 when second, and so on
    repeatTrials        = find(repeatNum>1);                    % gets indices of repeat trials
    
    conditions          = diff(b.stimuli');   % - when high freq & vis stim on left >> need to turn right, + when low freq & vis stim right and need to turn left
    conditions          = conditions(1:ntr);
    stimLeft            = find(conditions<0);
    stimRight           = find(conditions>0);
    
    
    firstAttempts               = response; % -1 when turned right (ie like high stim), 1 when turned left (ie like low stim)
    firstAttempts(repeatTrials) = NaN;
    firstAttempts               = firstAttempts(isfinite(firstAttempts));
    % now put it into format that can go into D for GLM model
    r = firstAttempts';
    r(firstAttempts==0) = 3;
    r(firstAttempts==1) = 2;
    r(firstAttempts==-1)= 1;
    
    cond                = conditions;
    cond(repeatTrials)  = NaN;
    cond                = cond(isfinite(cond));
    % put it into format that it can go into D for GLM model
    c = zeros(length(cond),2);
    cl = find(cond<0); cr = find(cond>0);
    c(cl,1) = abs(cond(cl));
    c(cr,2) = cond(cr);
    
    f = [b.evts.feedbackValues(1:ntr)]';
    f(repeatTrials) = NaN;
    f = f(isfinite(f));
    
    rN = [b.evts.repeatNum(1:ntr)]';
    rN(repeatTrials) = NaN; rN = rN(isfinite(rN));
    rN(rN==1) = 0;
    rN(rN>1) = 1;
    
    vv = vv+1;
    if isempty(D.response)
        % for GLM model
        D.contrast_cond = c;
        D.response      = r;
        D.feedbackType  = f;
        D.repeatNum     = rN;
        % to get errorbars using EJ's function
        stim_cond_EJ    = cond;
        response_EJ     = firstAttempts';
    else
        % for GLM model
        D.contrast_cond = cat(1,D.contrast_cond,c);
        D.response      = cat(1,D.response,r);
        D.feedbackType  = cat(1,D.feedbackType,f);
        D.repeatNum     = cat(1,D.repeatNum,rN);
        % to get errorbars using EJ's function
        stim_cond_EJ    = cat(2,stim_cond_EJ,cond);
        response_EJ     = cat(1,response_EJ,firstAttempts');
        
    end
end

%% make figure
unique_conditions  = unique(stim_cond_EJ);
cLeft = find(unique_conditions<0);
cZero = find(unique_conditions==0);
cRight= find(unique_conditions>0);
ucl = length(unique_conditions);

[~,~,~,percent_right,percent_left,percent_NoGo,ste_right,ste_left,ste_NoGo,~,~,~] = ...
    psych_curve_nogo(response_EJ,stim_cond_EJ,unique_conditions,ucl);

gv = GLM(D);
gv = gv.setModel('C50-subset');
gv = gv.fit;

maxC = max(max(gv.data.stimulus));
evalC = [linspace(maxC,0,100)', zeros(100,1);
    zeros(100,1), linspace(0,maxC,100)'];
evalC1d = evalC(:,2) - evalC(:,1);
phat = gv.calculatePhat(gv.parameterFits,evalC);
phatOrdered = zeros(size(phat));
phatOrdered(:,1) = phat(:,2);
phatOrdered(:,2) = phat(:,3);
phatOrdered(:,3) = phat(:,1);

spInds = [3 2 1];
for sp = 1:3
    subplot(1,3,spInds(sp))
    
    plot([0 0],[0 1],'k:'); hold on;
    plot([-1 1],[0.5 0.5],'k:');
    
    switch sp
        case 1
            errorbar(unique_conditions(cRight), percent_left(cRight), ste_left(cRight),...
                'ko', 'markerfacec',[0 0.9 0.4],'markersize',12);
            errorbar(unique_conditions(cZero), percent_left(cZero), ste_left(cZero),...
                'ko', 'markerfacec',[0.64 0.08 0.18],'markersize',12);
            errorbar(unique_conditions(cLeft), percent_left(cLeft), ste_left(cLeft),...
                'ko', 'markerfacec',[0.64 0.08 0.18],'markersize',12);
        case 3
            errorbar(unique_conditions(cRight), percent_right(cRight), ste_right(cRight),...
                'ko', 'markerfacec',[0.64 0.08 0.18],'markersize',12);
            errorbar(unique_conditions(cZero), percent_right(cZero), ste_right(cZero),...
                'ko', 'markerfacec',[0.64 0.08 0.18],'markersize',12);
            errorbar(unique_conditions(cLeft), percent_right(cLeft), ste_right(cLeft),...
                'ko', 'markerfacec',[0 0.9 0.4],'markersize',12);
        case 2
            errorbar(unique_conditions(cRight), percent_NoGo(cRight), ste_right(cRight),...
                'ko', 'markerfacec',[0.3 0.3 0.3],'markersize',12);
            errorbar(unique_conditions(cZero), percent_NoGo(cZero), ste_right(cZero),...
                'ko', 'markerfacec',[0 0.9 0.4],'markersize',12);
            errorbar(unique_conditions(cLeft), percent_NoGo(cLeft), ste_right(cLeft),...
                'ko', 'markerfacec',[0.3 0.3 0.3],'markersize',12);
    end
    
    %     gv.plotFit; hold on;
    plot(evalC1d,phatOrdered(:,sp),'Color',[0.5 0.5 0.5],'LineWidth',1);
    
    xlim([-0.055 0.055]);
    ylim([0 1]);
    box off
    xticks([-0.03 0.03])
    xticklabels({'High Tones','Low Tones'});
    xlabel('Contrast (%)');
    switch sp
        case 1
            ylabel('P(Low Tone Choice)');
        case 2
            ylabel('P(NoGo)');
        case 3
            ylabel('P(High Tone Choice)');
    end
    set(gca, 'FontSize', 12);
    axis square;
    
end

x0=1200; % 200 700 1200
y0=200;
width = 350;
height= 350;
set(gca,'units','points','position',[x0,y0,width,height]);
% xlim([-0.55 0.55]);
box off
set(gca, 'FontSize', 24)

