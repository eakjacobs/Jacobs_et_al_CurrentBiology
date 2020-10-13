function [rGO, rNOGO, rGoCorr, rGoInc, rCNG, rGOc, rInZC] = getChoiceNeglectInds_separateZC(b,ntr,stimName)
% function that returns indices of 
% GO     = choice, correct and incorrect
% NOGO   = incorrect nogo)
% GoCorr = correct choices
% GoInc  = incorrect choices
% CNG    = correct nogo
% rGOc   = GO excluding incorrect zero contrasts
% InZC    = incorrect zero contrast
% trials within behavioural block b
%
% written by Elina Jacobs, UCL Cortexlab

if nargin<3
    if ~isfield(b.params,'stimType')        % for concatenated blocks
        stimName = b.params{1,1}.stimType{1};
    else
        stimName = b.params.stimType{1};
    end
end

ee = [b.evts];

if isfield(b,'excludeFirstTrial');
    if b.excludeFirstTrial
        rV = [ee.responseValues(2:ntr+1)];
        fV = [ee.feedbackValues(2:ntr+1)];
    else
        rV = [ee.responseValues(1:ntr)];
        fV = [ee.feedbackValues(1:ntr)];
    end
else
    rV = [ee.responseValues(1:ntr)];
    fV = [ee.feedbackValues(1:ntr)];
    warning('if this is one of the early datasets, this is most likely going to be faulty!');
end

% get stimuli for viewing stim responses
if isfield(b.params,'stimType')
    stimType = b.params.stimType{2};
else
    stimType = b.params{1}.stimType{2};
end
if isfield(b.stimuli,'visContrasts')    % for audiovisual blocks
    stimuli = b.stimuli.visContrasts(1:ntr,:);
else
    stimuli = b.stimuli(1:ntr,:);
end
if isfield(b,'excludeFirstTrial');
    if b.excludeFirstTrial
        stimuli = b.stimuli(2:ntr+1,:);
    end
end
% find zero contrast trials
c0L = find(stimuli(:,1)==0);
c0R = find(stimuli(:,2)==0);
c0  = intersect(c0L,c0R);

rGO  = find(rV~=0);

rNoGo = find(rV==0);
fNoGo = fV(rNoGo);          % gives what feedback animal received for giving NoGo response
fNoGoPos = find(fNoGo>0);
fNoGoNeg = find(fNoGo<0);
rNOGO = rNoGo(fNoGoNeg);
rCNG  = rNoGo(fNoGoPos);

r0 = rV(c0);
r0Neg = find(r0~=0);
rInZC = c0(r0Neg);
% this way, incorrect trials in response to zero contrast are included in this vector as well as the choice vector

rGOc  = setdiff(rGO,rInZC);  % now, incorrect zero contras trials excluded from choice trials


switch stimName
    case 'Auditory'
        switch b.animal
            case 'EJ015'
                rT = b.evts.responseTimes - b.evts.interactiveOnT;
                zz = find(rT>7.5);
                rNOGO=zz;
                rGO(zz) = [];
        end
end

fGo  = fV(rGO);
fGoPos = find(fGo>0);
fGoNeg = find(fGo<0);

rGoCorr = rGO(fGoPos);
rGoInc  = rGO(fGoNeg);
rGoInc  = setdiff(rGoInc,rInZC);

end