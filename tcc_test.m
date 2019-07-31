%% Generate basis set and properties

% Create a bunch of basis functions
xs = 0:pi/64:2*pi;
basisSetMeans = pi/64:pi/64:2*pi;
k = exp(1);
% Plot the basis function responses over the 0:2pi range
figure(1); clf; hold on;
for i = 1:length(basisSetMeans)
%     plot(xs,exp(-k*abs(xs-basisSetMeans(i))));
    plot(xs,vonMises(xs,basisSetMeans(i),2));
end

%% Using the basis set, compute the correlation between noisy channel responses and ideal responses
% check if this reconstructs the original data?

clear aCR
for xi = 1:length(xs)
    aCR(:,xi) = vonMises(xs(xi),basisSetMeans,1);
end

clear thetas

for i = 1:10000
    % generate some data for theta = pi
    theta = pi;
    idealChannelResponses = vonMises(theta,basisSetMeans,1);
    idealChannelResponses = idealChannelResponses ./ max(idealChannelResponses);
    noisyChannelResponses = 1*idealChannelResponses + randn(size(idealChannelResponses));
    plot(noisyChannelResponses);
    ccorr = noisyChannelResponses*aCR; % just dot product
%     ccorr = corr(noisyChannelResponses',aCR);
    midx = find(ccorr==max(ccorr),1);
    thetas(i) = xs(midx);
end

figure(1); clf; hold on
[b,x] = hist(thetas,xs);

b = b ./ sum(b);
plot(x,b,'ok');

lx = 0:pi/64:pi;
like = computeTCCPDF(lx,1);
lx = [-fliplr(lx) lx];
like = [fliplr(like) like];
like = like ./ sum(like);
plot(lx+theta,like,'-r');
% plot(xs,

% figure; hold on
% plot(basisSetMeans,idealChannelResponses,'-b');
% plot(basisSetMeans,noisyChannelResponses,'-r');
% plot(xs,ccorr);
% 
% circcorr_xy = ifft(fft(idealChannelResponses').*conj(fft(noisyChannelResponses')));

%% Generate dissimilarity (1-corr)

% Compute the correlation between each basis set and the next one to
% compute the dissimilarity
figure;
startX = pi;
v = exp(-k*abs(startX-basisSetMeans));
dissimilarity = zeros(size(basisSetMeans));
rng = pi/64:pi/32:2*pi;
for ii = 1:length(v)
    v2 = circshift(v,[0,ii]);
    dissimilarity(ii) = 1-corr(v',v2');
end

% Normalize
dissimilarity = dissimilarity./max(dissimilarity);
%% Compare to Schurgin et al. data

data = [0, 0
6.998601810890314, 0.07794648395357928
15.961587671370282, 0.23378622650663594
24.382166230969382, 0.40601323673613965
36.75951761452152, 0.5741869315378648
50.404451277306144, 0.6931947290522409
69.267844311745, 0.7917701319784072
91.36615154077165, 0.9006282638913059
111.67510425004454, 0.9602856973237885
134.50743005471978, 0.9871849724354083
157.52215505364504, 0.9997420617455053
174.78154062834577, 1
179.6359385779414, 1];


figure(2); clf; hold on
plot(data(:,1),data(:,2),'o');
plot(rng*180/pi,dissimilarity,'-k');
axis([0 180 0 1]);

%% Generate choice data from exponential functions
% If the above is true, then the readout is operating over a bunch of
% exponentially distributed channels. With exponent ~ equal to exp(1) 
figure
% create the true basis set
n = 100;
basisSetMeans = pi/n:pi/(n/2):2*pi;

basis = repmat(xs,n,1);
basis = basis - repmat(basisSetMeans',1,length(xs));
basis = exp(-k*abs(basis));

imagesc(basis);





%% Explore the effects of attention in a channel model

% first generate a set of basis channels

% range over which to evaluate functions
xs = 0:pi/64:2*pi;
% basis set
basisSetMeans = pi/64:pi/8:2*pi;

% each channel has a von mises tuning function, we will now simulate 1000
% repeats of generating channel responses and "reading out" the true
% direction (theta=pi). 

% aCR holds all the channel responses 
clear aCR
for xi = 1:length(xs)
    aCR(:,xi) = vonMises(xs(xi),tempMeans,1);
end

clear thetas

% temporarily move basisSetMeans closer together
tempMeans = basisSetMeans - (basisSetMeans-pi)/3;

for i = 1:10000
    % generate some data for theta = pi
    theta = pi;
    idealChannelResponses = vonMises(theta,tempMeans,1);
    idealChannelResponses = idealChannelResponses ./ max(idealChannelResponses);
    noisyChannelResponses = idealChannelResponses + randn(size(idealChannelResponses));
    plot(noisyChannelResponses);
    ccorr = noisyChannelResponses*aCR; % just dot product
%     ccorr = corr(noisyChannelResponses',aCR);
    midx = find(ccorr==max(ccorr),1);
    thetas(i) = xs(midx);
end

% this generates a bunch of thetas with some distribution. 




% same thing but applying a gain to some channels 

%% 
xs = 0:pi/16:2*pi;
cmap = colorblindmap/255;
% using nthetas, athetas, cthetas
h = figure; hold on

[cb,~] = hist(cthetas,xs);

plot(xs,cb,'-','Color',cmap(7,:));

[ab,~] = hist(athetas,xs);

plot(xs,ab,'-','Color',cmap(6,:));

[nb,~] = hist(nthetas,xs);

plot(xs,nb,'-k');

set(gca,'XTick',[0 pi 2*pi],'XTickLabel',[-360 0 360]);
set(gca,'YTick',[0 1000]);

xlabel('Distance from target (deg)');
ylabel('Response density (pdf)');

legend({'Attended: Shift','Attended: Gain','Unattended'});
drawPublishAxis('figSize=[4.5,4.5]');

savepdf(h,fullfile('~/proj/afcom/figures/','TCC_channel_attention_likelihoods.pdf'));










