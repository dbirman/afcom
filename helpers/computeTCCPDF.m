function like = computeTCCPDF(rads,dprime)
%% COMPUTE TCCPDF
% Usage: likelihood = computeTCCPDF(xs,dprime)
%
% Inputs:
%   rads - range of x values to evaluate at in *radians*
%   dprime - d prime value of the model
%
% Outputs:
%   likelihood - the likelihood that the TCC model will sample each x,
%                normalized

N = 100; % number of channels, anything over 100 converges 

if any(rads>pi)
    disp('You might have submitted degrees, converting to radians');
    rads = rads * pi / 180;
end

if dprime==0
    like = ones(size(rads));
    like = like./sum(like);
    return
end

if iscolumn(rads)
    rads = rads';
end
% set up the encoders
px = 1-pscale(abs(rads)*180/pi);
% px = px * dprime;
sigma = 1 / dprime;
% set up the encoder response-range
minmax = 4*sigma;
rrange = linspace(-minmax,1+minmax,N);
dr = rrange(2)-rrange(1);

% pre-compute the probability density functions
like = zeros(1,length(px));

pxidxs = 1:length(px);

% setup the normpdf in advance
allpdfs = normpdf(repmat(rrange,length(px),1),repmat(px',1,length(rrange)),sigma)*dr;

for xi = pxidxs
    % for each x location, compute the probability that the encoder at this
    % location exceeds all of the others
    nx1 = px([1:(xi-1) (xi+1):end]);    
    like(xi) = allpdfs(xi,:) * prod(normcdf(repmat(rrange',1,length(nx1)),repmat(nx1,length(rrange),1),sigma),2);
end

if sum(like)<0.99
    like = eps*ones(size(like));
    warning('Sum of TCC likelihood is not close enough to 1, the range may not be appropriate');
elseif sum(like>1.01)
    warning('Error');
    stop = 1;
end