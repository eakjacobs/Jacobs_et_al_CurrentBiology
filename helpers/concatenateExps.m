function [allExps] = ...
    concatenateExps(mousename, expDate, expList, doHemoCorrect, fewerSVDs, useDFF, noHPfilter)
% function that concatenates all behavioural experiments from one day

% set some defaults
if nargin < 7
    noHPfilter    = false;
end
if nargin < 6
    useDFF      = false;
end
if nargin < 5
    fewerSVDs     = false;
end
if nargin < 4
    doHemoCorrect = true;
end

SetDefaultDirs;
BehDir      = DIRS.expInfo;
svdDirOld   = '\\zserver.cortexlab.net\Data\GCAMP';
svdDirNew   = fullfile(serverName,'Data','Subjects');

%% get experiments of that day, see which ones were recorded
if isempty(expList)
    thisDir = dir(fullfile(BehDir, mousename, expDate));
    thisDir = thisDir(3:end);
    
    expNums = str2double({thisDir.name});
    
    % check whether they were all recorded, set the ones that weren't to NaN
    for iie=1:length(expNums)
        thisExp = num2str(expNums(iie));
        Vnpy = fullfile(svdDirOld,mousename,expDate,thisExp,strcat(expDate,'_',thisExp,'_',mousename,'_SVD_V.npy'));
        svdDir = 'old';
        if ~exist(Vnpy)
            Vnpy = fullfile(svdDirOld,mousename,expDate,thisExp,'SVD_Results_V.npy');
            if ~exist(Vnpy)
                Vnpy = fullfile(svdDirNew,mousename,expDate,thisExp,'svdTemporalComponents_blue.npy');
                svdDir = 'new';
                if ~exist(Vnpy)
                    Vnpy = fullfile(svdDirNew,mousename,expDate,'svdTemporalComponents_blue.npy');
                    if ~exist(Vnpy)
                        expNums(iie) = NaN;
                    end
                end
            end
        end
    end
    clear iie
    expList = expNums(isfinite(expNums));
end

allB = [];
allB.animal = mousename;
allB.iseries = expDate;
allB.exps   = expList;
allB.trialsPerExp = [];
allB.completedTrials = [];

allV = [];
allU = [];
allT = [];

nv = 0;
na = 0;
nav = 0;

for iie = 1:length(expList)
    thisExp = num2str(expList(iie));
    
    %% get behavioural data
    expRef  = strcat(expDate,'_',thisExp,'_',mousename);
    Exps.animal = mousename;
    Exps.iseries = expDate;
    Exps.iexp   = thisExp;
    [b] = generateGenBlock(expRef, Exps);
    
    ntr = b.completedTrials;
    allB.trialsPerExp(iie) = ntr;
    if isempty(allB.completedTrials)
        allB.completedTrials = ntr;
    else allB.completedTrials = allB.completedTrials + ntr;
    end
    
    if ~isfield(allB,'evts')
        evFields = fieldnames(b.evts);
        allB.evts = cell2struct(cell(size(evFields)), evFields);
        
        for iif = 1:length(evFields)
            allB.evts.(evFields{iif}) = b.evts.(evFields{iif})(1:ntr);
        end
        
        allB.wheel = b.wheel;
    else
        for iif = 1:length(evFields)
            switch evFields{iif}(end-4:end)
                case 'Times'
                    b.evts.(evFields{iif})    = bsxfun(@plus,b.evts.(evFields{iif}),prevExpEndTime);
            end
            switch evFields{iif}(end)
                case 'T'
                    b.evts.(evFields{iif})    = bsxfun(@plus,b.evts.(evFields{iif}),prevExpEndTime);
            end
            allB.evts.(evFields{iif}) = cat(2,allB.evts.(evFields{iif}),b.evts.(evFields{iif})(1:ntr));
        end
        
        allB.wheel.Values = cat(2,allB.wheel.Values,b.wheel.Values);
        allB.wheel.Times  = cat(2,allB.wheel.Times,b.wheel.Times+prevExpEndTime);
    end
    
    if ~isfield(allB,'stimuli')
        if isfield(b.stimuli,'visContrasts')        % for AV block case
            allB.stimuli = b.stimuli.visContrasts(1:ntr,:);
        else
            allB.stimuli = b.stimuli(1:ntr,:);
        end
    else
        if isfield(b.stimuli,'visContrasts')        % for AV block case
            allB.stimuli = cat(1,allB.stimuli,b.stimuli.visContrasts(1:ntr,:));
        else
            allB.stimuli = cat(1,allB.stimuli,b.stimuli(1:ntr,:));
        end
    end
    
    allB.params{iie} = b.params;
    allB.contRecording = b.contRecording;
    allB.expType = b.expType;
    allB.rigName = b.rigName;
    
    
    %% get neural data
    if noHPfilter
        [U, V, t, SVDinfo]  = loadSVDfiles_withoughtHighPassFilter(b, doHemoCorrect, fewerSVDs);
    else
        [U, V, t, SVDinfo]  = loadSVDfiles(b, doHemoCorrect, fewerSVDs);
    end
    
    if useDFF
        disp('applying DFF...');
        [U, V] = applyDFF(U,V,b);
    end
    
    if isempty(allT);
        allT = t;
        allV = V;
    else
        allT = cat(2,allT,t+prevExpEndTime);
        allV = cat(2,allV,V);
    end
    switch SVDinfo
        case 'GCAMPraw'
            for iff = 1:length(expList)
                Ufn{iff} = strcat('exp',num2str(iff));
            end
            allU.(Ufn{iie}) = U;
    end
    
    prevExpEndTime = allT(end)+10;
    tOffset(iie) = prevExpEndTime;
    
    allSVDinfo{iie} = SVDinfo;
    
    
end

% given that there is only one U for the newer exps, this should work
if isempty(allU)
    allU = U;
end

allExps.block = allB;
allExps.U     = allU;
allExps.V     = allV;
allExps.t     = allT;
allExps.tOffset = tOffset(1:end-1);
allExps.SVDinfo = allSVDinfo;

end