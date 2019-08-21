%%

xs = 0:pi/64:2*pi;

basisSetMeans = 0:pi/256:2*pi;

for bi = 1:length(basisSetMeans)
    aCR(bi,:) = vonMises(xs,basisSetMeans(bi),20);
end

figure;
plot(aCR');

%% 
dprime = 2;

% compute a sample

theta_est = [];

for i = 1:10000
    theta = pi;
    sample = dprime*vonMises(basisSetMeans-theta,0,1) + randn(1,length(basisSetMeans));
    
    % least squares computation
%     dist = sum((aCR - repmat(sample',1,size(aCR,2))).^2);
%     theta_est(i) = xs(find(dist==min(dist),1));
%     
    % dot product computation
    dist = sample * aCR;
    theta_est(i) = xs(find(dist==max(dist),1));
end

figure(1); clf; hold on

tcc_like = computeTCCPDF(0:pi/64:pi,dprime);
tcc_like = [fliplr(tcc_like) tcc_like(2:end)];
tcc_like = tcc_like./sum(tcc_like);

[b,~] = hist(theta_est,xs);
b = b./sum(b);

plot([-pi:pi/64:pi]+pi,tcc_like,'-k');
plot(xs,b,'-r');

temp = vonMises(xs-pi,0,dprime*3);
temp = temp ./ sum(temp);
plot(xs,temp);