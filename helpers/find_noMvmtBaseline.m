function [baseline] = find_noMvmtBaseline(baseline, ntr, Fs);
% finds last moment of movement during the baseline (allowing for a little
% bit of wiggling) prior to stimulus/cue onset
%
% then adds field to baseline structure that contains the frames per trial
% in which no movement occurred
%
% written by Elina Jacobs, UCL Cortexlab

baseline.noMvmtFrames       = zeros(ntr,2);
baseline.noMvmtStartIndex   = zeros(ntr,1);     
% gives index within each trial when the movement stops within
% wheelValues/pupilTraces

baseline.noMvmtFrames(:,2)  = baseline.V_OnOffFramenumbers(:,2);
nf = size(baseline.wheelValues,2);
ww = baseline.wheelValues;
xx = find(ww>0&ww<=1|ww<0&ww>=-1);        % allow for tiny twitches at the end
ww(xx) = 0;

for iiw=1:ntr
    lm = find(ww(iiw,:),1,'last');
    if isempty(lm)
        lm = 0;
    end
    baseline.noMvmtFrames(iiw,1) = baseline.noMvmtFrames(iiw,2) - (nf-lm);
    % if mouse twitched at stim onset, that would make only the very last value zero
    if diff(baseline.noMvmtFrames(iiw,:)') == 1       
        thisW = ww(iiw,:);
        wv = thisW(end-1);
        thisW = thisW-wv;
        thisW(end) = 0;
        xx2 = find(thisW>0&thisW<=1|thisW<0&thisW>=-1);
        thisW(xx2) = 0;
        lm = find(thisW,1,'last');
        baseline.noMvmtFrames(iiw,1) = baseline.noMvmtFrames(iiw,2) - (nf-lm);
    end
    baseline.noMvmtStartIndex(iiw) = lm+1;
    % if baseline with no movement too short, allow for a bit more
    % twitching
    if diff(baseline.noMvmtFrames(iiw,:)') < 0.7*Fs
        if diff(baseline.noMvmtFrames(iiw,:)') == 1
            thisW = ww(iiw,:);
            wv = thisW(end-1);
            thisW = thisW-wv;
            thisW(end) = 0;
            xx2 = find(thisW>0&thisW<=1|thisW<0&thisW>=-1);
            thisW(xx2) = 0;
        else
            thisW = baseline.wheelValues(iiw,:);
        end
        thisW = fliplr(thisW(isfinite(thisW)));
        dw  = diff(thisW);
        xx2 = find(dw>0&dw<=2|dw<0&dw>=-2);
        dw(xx2) = 0;
        lm2 = find(dw,1,'first');
        if isempty(lm2)
            lm2 = length(thisW);
        end
        baseline.noMvmtFrames(iiw,1) = baseline.noMvmtFrames(iiw,2) - lm2;
        baseline.noMvmtStartIndex(iiw) = nf - lm2;
    end
    % if baseline with no movement still too short, again allow for a bit more
    % twitching
    if diff(baseline.noMvmtFrames(iiw,:)') < 0.7*Fs
        if diff(baseline.noMvmtFrames(iiw,:)') == 1
            thisW = ww(iiw,:);
            wv = thisW(end-1);
            thisW = thisW-wv;
            thisW(end) = 0;
            xx2 = find(thisW>0&thisW<=1|thisW<0&thisW>=-1);
            thisW(xx2) = 0;
        else
            thisW = baseline.wheelValues(iiw,:);
        end
        thisW = fliplr(thisW(isfinite(thisW)));
        dw  = diff(thisW);
        xx2 = find(dw>0&dw<=3.5|dw<0&dw>=-3.5);
        dw(xx2) = 0;
        lm2 = find(dw,1,'first');
        if isempty(lm2)
            lm2 = length(thisW);
        end
        baseline.noMvmtFrames(iiw,1) = baseline.noMvmtFrames(iiw,2) - lm2;
        baseline.noMvmtStartIndex(iiw) = nf - lm2;
    end
end

baseline.noMvmtStartIndex(baseline.noMvmtStartIndex==0) = 1;

wd = diff(baseline.noMvmtFrames');
if min(wd) < 0.7*Fs
    warning('there are some no movement periods that are too short for a meaningful analysis');
end

end