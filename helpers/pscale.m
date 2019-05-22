function y = pscale(x)

if all(abs(x)<pi)
    warning('You called p-scale and passed in what appear to be radians, converting to degrees');
    x = x*180/pi;
end
% NAKA RUSHTON VERSION

rmax = 1.1;
n = 1.5; 
c50 = 35;
b = 0;

y = rmax * x.^n ./ (x.^n + c50^n) + b;

% POWER LAW VERSION

% k = .02;
% a = 180*k;
% 
% if ~inv
% 
%     % x = 0:180;
%     y = pscale_(x,a,k);
%     y = y./pscale_(180,a,k);
%     
% else
%     
%     y = (log(x-a) + log(a))/-k;
% end

return

function y = pscale_(x,a,k)

y = a-(a*exp(-k*x));

return

%% data from paper figure
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

%% power law test
figure(1);

k = .03;
a = 180*k;

x = 0:180;

y = a-(a*exp(-k*x));
y = y./max(y);

clf
hold on
plot(data(:,1),data(:,2),'or');
plot(x,y);
plot([0 180],[0 1],'--k');

%% naka-rushton version
rmax = 1.1;
n = 1.5; 
c50 = 35;
b = 0;

y = rmax * x.^n ./ (x.^n + c50^n) + b;
h = figure;
clf
hold on
plot(x,y,'-k');
plot(data(:,1),data(:,2),'o','MarkerFaceColor','k','MarkerEdgeColor','w');
xlabel('Distance between reference and mean of pair (deg)');
ylabel('Psychological distance');
axis([0 180 0 1]);
set(gca,'XTick',[0 90 180]);
set(gca,'YTick',[0 1]);
legend({'Fit','Data from Schurgin et al.'});
drawPublishAxis('figSize=[16,9]','poster=1');
% plot([0 180],[0 1],'--k');
savepdf(h,fullfile('~/proj/afcom/figures/psychdist.pdf'));

%% test exponential version
k = 40;
y = 1-exp(-x/k);

h = figure(4);
clf
hold on
plot(x,y,'-k');
plot(data(:,1),data(:,2),'o','MarkerFaceColor','k','MarkerEdgeColor','w');
xlabel('Distance between reference and mean of pair (deg)');
ylabel('Psychological distance');
axis([0 180 0 1]);
set(gca,'XTick',[0 90 180]);
set(gca,'YTick',[0 1]);
legend({'Fit','Data from Schurgin et al.'});
% drawPublishAxis('figSize=[16,9]','poster=1');
% plot([0 180],[0 1],'--k');
% savepdf(h,fullfile('~/proj/afcom/figures/psychdist.pdf'));