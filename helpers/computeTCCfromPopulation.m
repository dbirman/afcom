function like = computeTCCfromPopulation(rads,sigma,p)
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

if iscolumn(rads)
    rads = rads';
end



npop = 100;
rpeaks = ((2*pi)/npop/2-pi):(2*pi/npop):pi;

responses = zeros(length(rads),length(rpeaks));
for ri = 1:length(rads)
    responses(ri,:) = normpdf(rpeaks,rads(ri),sigma);
end

% weight the responses by the exponential
weights = exp(-p*abs(rpeaks));
weights = weights ./ sum(weights);

like = responses * weights';

stop = 1;

%% test script
% rads = -pi:pi/128:pi;
rads = alldata(:,4);

% pull out histogram
xs = 0:pi/32:pi;
[n,xs] = hist(rads,xs);

p = 0.1;
sigma = pi/8;


npop = 100;
rpeaks = ((2*pi)/npop/2-pi):(2*pi/npop):pi;

responses = zeros(length(xs),length(rpeaks));
for ri = 1:length(xs)
    responses(ri,:) = normpdf(rpeaks,xs(ri),sigma);
end

% weight the responses by the exponential
weights = normpdf(rpeaks,0,p);%rpeaks./sum(rpeaks);%exp(-p*abs(rpeaks));
% weights = weights ./ sum(weights);

like = responses * weights';

figure(1);
clf
hold on

plot(xs,n./sum(n),'ok');
% plot(rpeaks,response,'ok');
plot(xs,like./sum(like),'-r');

% now plot the TCC model prediction
plot(xs,computeTCCPDF(xs,1.6),'-b');

legend({'Data','Population readout','TCC'});