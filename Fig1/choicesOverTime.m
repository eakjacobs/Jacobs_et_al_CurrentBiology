function [percCorrect, percChoice, percIncorrect, percNogo, percIncorrNogo] = choicesOverTime(b,ww);
% computes percent correct, incorrect, choice (which is correct and
% incorrect choices), nogo and incorrect nogo over time within one
% behavioural session
%
% inputs:
% b - block (in format from generateGenBlock)
% ww - window width (over how many trials to average)

if nargin<2
    ww = 20;        % defines window width
end

[EJDirs] = setEJDirs;

ntr = b.completedTrials;
ee = [b.evts];

if b.excludeFirstTrial
    rV = [ee.responseValues(2:ntr)];
    fV = [ee.feedbackValues(2:ntr)];
    ntr = ntr-1;
else
    rV = [ee.responseValues(1:ntr)];
    fV = [ee.feedbackValues(1:ntr)];
end

%% first find percent correct over time

fVc = b.evts.feedbackValues;
fVc(fVc<0) = 0;

percCorrect = performanceOverTime(ww,ntr,fVc);


%% find percent choices over time

rGo  = find(rV~=0);
percGo = zeros(1,ntr);
percGo(rGo) = 1;

percChoice = performanceOverTime(ww,ntr,percGo);

%% find percent nogo over time, divide into correct and incorrect as appropriate

rNoGo = find(rV==0);
fNoGo = fV(rNoGo);          % gives what feedback animal received for giving NoGo response
fNoGoPos = find(fNoGo>0);
fNoGoNeg = find(fNoGo<0);

if isempty(fNoGoPos)        % giving a NoGo response is always incorrect

    pInogo = zeros(1,ntr);
    pInogo(rNoGo) = 1;
    percIncorrNogo = performanceOverTime(ww,ntr,pInogo);
    
    percNogo = percIncorrNogo;
    
else
    
    pNogo = zeros(1,ntr);
    pNogo(rNoGo) = 1;
    percNogo = performanceOverTime(ww,ntr,pNogo);
    
    pInogo = zeros(1,ntr);
    pInogo(rNoGo(fNoGoNeg)) = 1;
    percIncorrNogo = performanceOverTime(ww,ntr,pInogo);
    
end

%% find percent incorrect choices over time

fVi = b.evts.feedbackValues;
fVi(fVi>0) = 0;
fVi(rNoGo) = 0;
fVi = abs(fVi);

percIncorrect = performanceOverTime(ww,ntr,fVi);

end

