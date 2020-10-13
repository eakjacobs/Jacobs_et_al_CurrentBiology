function [dffU, dffV] = applyDFF(U,V,b);
% function that loads meanImage & then calls dffFromSVD function to
% generate dffU & dffV

% written by Elina Jacobs, UCL Cortexlab

SetDefaultDirs;
svdDir= fullfile(serverName,'Data','Subjects');
svdDirOld   = '\\zserver.cortexlab.net\Data\GCAMP';
svdDirNew   = '\\zserver.cortexlab.net\Data\Subjects';

if exist(fullfile(svdDirNew,b.animal,b.iseries,b.iexp),'dir')
    avgIm = readNPY(fullfile(svdDirNew,b.animal,b.iseries,'meanImage_blue.npy'));
else
    load(fullfile(svdDirOld,b.animal,b.iseries,b.iexp,'dataSummary.mat'));
    avgIm = dataSummary.meanImage(ops.yrange,ops.xrange);
end

if b.animal(1:2) == 'EJ'
    % rotate U & meanImage so that visual cortex at bottom, auditory on left, somsen on top
    if str2num(b.iseries(4)) == 5   % if the experiment was done in 2015, need to do the following to get visual cortex to bottom, auditory to left etc
        avgIm = rot90(rot90(avgIm));
    else
        if str2num(b.iseries(4)) > 5        % if the experiment was done after 2015
            avgIm      = flipud(avgIm);
            if str2num(b.iseries(6:7)) < 8 % if the experiment was done before August 2016
                avgIm      = rot90(avgIm);
            end
        end
    end
end

% if there is a small discrepency between average image and U, just cut
% avgIm to the size of U
if size(avgIm,1) > size(U,1)
    if size(avgIm,1)-size(U,1) < 3
        avgIm = avgIm(1:size(U,1),:);
    end
end
if size(avgIm,2) > size(U,2)
    if size(avgIm,2)-size(U,2) < 3
        avgIm = avgIm(:,1:size(U,2));
    end
end

[dffU, dffV] = dffFromSVD(U, V, avgIm);         
% dffFromSVD is a function from Cortexlab Widefield GitHub repository

end