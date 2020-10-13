function P = pixSpec_per_cond(V, U, FS, frb)
% Inputs:
% V         = matrix of nSVD x timepoints
% U         = matrix of xPixels x yPixels x nSVD
% FS        = sampling Frequency
% frb       = frequency band of interest: [lowFreq highFreq]
% Outputs:
% P = the power of each pixel in the frequency band 
%
% written by Michael Okun, UCL Cortexlab
% requires Chronux toolbox

assert(numel(size(U)) == 3 & size(U,3) == size(V, 1), 'inputs of wrong format')

V = bsxfun(@minus, V, mean(V, 2)); % mean-subtract each SVD

x = size(U,1);
y = size(U,2);
U = reshape(U, x*y, []); % pixels x nSVD

N = size(V,2);
nfft = max(2^nextpow2(N),N);
tapers = dpsschk([1 1], size(V,2),FS); % 1 taper
J = squeeze(mtfftc(V', tapers, nfft, FS)); % Fourier transform (frequencies x nSVD)
[~,findx] = getfgrid(FS,nfft,frb); 
J = J(findx, :); % leave only the frequencies we're interested in...

F = J*U'; % frequencies x pixels
P = sum(conj(F).*F)'; % the power 
P = reshape(P, x, y);

end