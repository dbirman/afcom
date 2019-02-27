%% Model of Schurgin et al. paper 

% assume five colors are encoded in memory randomly
clear aenc
clear choices
for rep = 1:1000
    
    cols = 180;

    dprime = 1;
    noise = 1;

    % for each color, 360 units encode it -- but according to the distance
    % function * the dprime

    nEncoders = 50;

    x = 0:360/nEncoders:360;

    % generate the distance matrix
    dist = abs(repmat(x',1,length(cols)) - repmat(cols,length(x),1));
    % rotate distances > 180
    dist(dist>180) = 180-mod(dist(dist>180),180);

    encoded = dprime * (1-pscale(dist,false)) + randn(length(x),length(cols))*noise;
    
    choice = find(encoded==max(encoded),1);
    choices(rep) = choice;
end

%% Plot

h = figure(1);
clf

vline(cols(1),'--r');

d = 3;

bins = (d/2):d:(360-d/2);

n = hist(x(choices),bins);
n = n./sum(n);
plot(bins,n,'o','MarkerFaceColor','k','MarkerEdgeColor','w');

% l = laplace(bins,180,30);
% l = l./sum(l);
% plot(bins,l,'-r');
% 
% vm = vonMises(bins*pi/180,pi,4);
% vm = vm./sum(vm);
% plot(bins,vm,'-b');

a = axis;
axis([0 360 0 0.05]);
% axis([0 360 0 20]);

sf = 3.5;

% trying to compute the probability distribution directly from bins
% like = zeros(size(bins));
% 
% for bi = 1:length(bins)
%     bin = bins(bi);
%     
% %     dist = abs(bin-cols);
%     dist = pscale(abs(bin-cols),false);
%     like(bi) = normpdf(dist,0,0.2);
% end
% 
% plot(bins,like);