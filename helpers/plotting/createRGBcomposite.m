function compI = createRGBcomposite(Im1,Im2,cAxLims,cm1)

% where Im1 is an image with colored colorscale,
% Im2 is an image with gray colorscale
% cAxLims gives the color axis limits; first row for Im1, second row for Im2
% cm1 gives colorscale (red white blue default)

% created by Elina Jacobs, UCL CortexLab

if nargin<4
    % get colourscale that is red for pos, blue for neg, and white for 0 values
    cm1 = RedWhiteBlue; cm1 = flipud(cm1);
end

if nargin<3
    cAxLims = [];
end


ncolors = size(cm1,1);                  % number of colours in the colorscale
cm2 = gray(ncolors);                    % make second black to white colorscale with same number of colors

Im1RGB = convertToRGB(Im1,cm1,ncolors,cAxLims(1,:));

Im2RGB = convertToRGB(Im2,cm2,ncolors,cAxLims(2,:));

compI  = Im1RGB.*Im2RGB;
