

%% example plot for poster
cmap = brewermap(15,'RdYlGn');
cmap = shift(cmap,[7,0]);
cmap = cmap(1:13,:);

% plot a TCC channel model, showing the activation of different units 
units = -pi:(2*pi/(size(cmap,1)-1)):pi;
units = sort(units);
x = -pi:pi/128:pi;

h = figure(2); clf; 
subplot(211); hold on
vline(0,'--k');
for ui = 1:length(units)
    mu = units(ui);
    profile = 1-pscale(abs(x-mu)*180/pi);
    % compute my activation at 0
    activation = 1-pscale(abs(mu)*180/pi);
    plot(x,profile,'-','Color',cmap(ui,:));
end
axis([-pi pi -0.1 1.1]);
set(gca,'XTick',[-pi 0 pi],'XTickLabel',[-180 0 180]);
set(gca,'YTick',[0 1]);
drawPublishAxis('figSize=[25,10]','poster=1');
subplot(212); hold on
for ui = 1:length(units)
    mu = units(ui);
    profile = 1-pscale(abs(x-mu)*180/pi);
    % compute my activation at 0
    activation = 1-pscale(abs(mu)*180/pi);
    errbar(mu,activation,1,'-','Color',cmap(ui,:));
    plot(mu,activation,'o','MarkerFaceColor',cmap(ui,:),'MarkerEdgeColor','w','MarkerSize',10);
end

ylabel('Channel activation (a.u.)');
xlabel('Distance from presented stimulus (deg)');
axis([-pi pi -1 2]);
set(gca,'XTick',[-pi 0 pi],'XTickLabel',[-180 0 180]);
set(gca,'YTick',[0 1]);
drawPublishAxis('figSize=[20,13]','poster=1');

savepdf(h,fullfile('~/proj/afcom/figures/','TCC model.pdf'));

%% example plot with von mises
h = figure;

subplot(121);
plot(x,repmat(0.1,size(x)),'-k');
axis([-pi pi 0 1]);
set(gca,'XTick',[-pi 0 pi],'XTickLabel',[-180 0 180]);
set(gca,'YTick',[0 1]);
drawPublishAxis('figSize=[20,5]','poster=1');
subplot(122);
plot(x,vonMises(x,0,2),'-k');
axis([-pi pi 0 1]);
set(gca,'XTick',[-pi 0 pi],'XTickLabel',[-180 0 180]);
set(gca,'YTick',[0 1]);
drawPublishAxis('figSize=[20,5]','poster=1');

savepdf(h,fullfile('~/proj/afcom/figures/','vonMises_example.pdf'));

%% example plot with different sensitivity (left panel scaled, right panel normalized)
x = -pi:pi/128:pi;

h = figure(2); clf; hold on

% cmap = colorblindmap/255;

cmap = brewermap(13,'Reds');

ax = pi*0.75;
ay = 0.4;
d = 0.5;

subplot(121); hold on
% compute the pscale function for four different angs
dprimes = logspace(-1,0.5,6);
dprimes = fliplr(dprimes);
legs = {};
for di = 1:length(dprimes)
    dprime = dprimes(di);
    
    activation = dprime*(1-pscale(abs(x)*180/pi));
    
    plot(x,activation,'-','Color',cmap(14-di,:));
    legs{end+1} = sprintf('d'' %1.2f',dprime);
end
legend(legs);


ylabel('Channel activation (sigma=1)');
xlabel('Distance from presented stimulus (deg)');
axis([-pi pi 0 3.2]);
set(gca,'XTick',[-pi 0 pi],'XTickLabel',[-180 0 180]);
set(gca,'YTick',[0 3.2]);
drawPublishAxis('figSize=[40,13]','poster=1');

subplot(122); hold on
% compute the pscale function for four different angs
dprimes = logspace(-1,0.5,6);
dprimes = fliplr(dprimes);
legs = {};
for di = 1:length(dprimes)
    dprime = dprimes(di);
    
    activation = computeTCCPDF(x,dprime);
    plot(x,activation,'-','Color',cmap(14-di,:));
    legs{end+1} = sprintf('d'' %1.2f',dprime);
end

ylabel('Response likelihood (pdf)');
xlabel('Distance from presented stimulus (deg)');
axis([-pi pi 0 0.08]);
set(gca,'XTick',[-pi 0 pi],'XTickLabel',[-180 0 180]);
set(gca,'YTick',[0 0.08]);
drawPublishAxis('figSize=[30,10]','poster=1');

savepdf(h,fullfile('~/proj/afcom/figures/','sensitivity_example.pdf'));

%%
x = -pi:pi/64:pi;
% also plot for these dprime values the tcc_data plots
for di = 1:length(dprimes)
    like = computeTCCPDF(x,dprimes(di));
    plot_tcc_data(x,like*4,cmap(14-di,:),sprintf('tcc_circle_%i.pdf',di));
end

%% TCC Channel model (for the afcom_avg task): showing the effect of all four directions of motion 
units = -pi:pi/128:pi;
units = sort(units);
x = -pi:pi/128:pi;

h = figure(2); clf; hold on

% cmap = colorblindmap/255;

cmap = [
               0 136 164
               0 136 164
               247 144 30
               247 144 30]/255;

types = {'target','side','feat','dist'};

ax = pi*0.75;
ay = 0.4;
d = 0.5;

% compute the pscale function for four different angs
angs = [0 -pi/2 -pi/4 pi/4];
weights = [0.5 0.5 0.2 0.2];
for ai = 1:length(angs)
    activation = 1-pscale(abs(angdist(units,angs(ai)))*180/pi);
    plot(x,weights(ai)*activation,'-','Color',cmap(ai,:));
    
    % compute the y scaling factor, which should be 
    ca = cos(pi/2-angs(ai));
    sa = sin(pi/2-angs(ai));
    yf = (1.3*sa)/(pi);
    
    ax1 = ax-d*ca;
    ax2 = ax+d*ca;
    ay1 = ay-d*sa*yf;
    ay2 = ay+d*sa*yf;
    arrow([ax1 ay1],[ax2,ay2],'Color',cmap(ai,:));
end
for ai = 1:length(angs)
    activation(ai,:) = 1-pscale(abs(units-angs(ai))*180/pi);
end
vline(mean(angs(1:2)),'--r');
plot(x,activation' * weights','-k');
ylabel('Channel activation (a.u.)');
xlabel('Distance from presented stimulus (deg)');
axis([-pi pi -0.3 1]);
set(gca,'XTick',[-pi 0 pi],'XTickLabel',[-180 0 180]);
set(gca,'YTick',[0 1]);
drawPublishAxis('figSize=[40,13]','poster=1');

savepdf(h,fullfile('~/proj/afcom/figures/','afcom_tcc_activationmodel_avg.pdf'));

%% TCC Channel model (for the afcom task): showing the effect of all four directions of motion 
units = -pi:pi/128:pi;
units = sort(units);
x = -pi:pi/128:pi;

h = figure(2); clf; hold on

cmap = ac_cmap;

types = {'target','side','feat','dist'};

ax = pi*0.75;
ay = 0.4;
d = 0.5;

% compute the pscale function for four different angs
angs = [pi rand*2*pi rand*2*pi rand*2*pi]-pi;
weights = [0.5 0.25 0.16 0.09];
for ai = 1:length(angs)
    activation = 1-pscale(abs(units-angs(ai))*180/pi);
    plot(x,weights(ai)*activation,'-','Color',cmap.(types{ai}));
    ca = cos(pi/2-angs(ai));
    sa = sin(pi/2-angs(ai));
    yf = (1.3*sa)/(pi);
    
%     ax1 = ax-d*ca;
%     ax2 = ax+d*ca;
%     ay1 = ay-d*sa*yf;
%     ay2 = ay+d*sa*yf;
%     arrow([ax1 ay1],[ax2,ay2],'Color',cmap.(types{ai}));
end
for ai = 1:length(angs)
    activation(ai,:) = 1-pscale(abs(units-angs(ai))*180/pi);
end
plot(x,activation' * weights','-k');
ylabel('Response likelihood (pdf)');
xlabel('Distance from presented stimulus (deg)');
axis([-pi pi 0 0.75]);
set(gca,'XTick',[-pi 0 pi],'XTickLabel',[-180 0 180]);
set(gca,'YTick',[0 0.75]);
drawPublishAxis('figSize=[25,10]','poster=1');

savepdf(h,fullfile('~/proj/afcom/figures/','afcom_tcc_activationmodel.pdf'));

%% Basic likelihood plot (as an example
units = -pi:pi/64:pi;
units = sort(units);
x = -pi:pi/64:pi;

h = figure(1); clf; hold on

% compute the pscale function for four different angs
like = computeTCCPDF(x,1);
plot(x,like,'-','Color','k');

ylabel('Response likelihood (pdf)');
xlabel('Distance from presented stimulus (deg)');
axis([-pi pi 0 0.05]);
set(gca,'XTick',[-pi 0 pi],'XTickLabel',[-180 0 180]);
set(gca,'YTick',[0 0.05]);
drawPublishAxis('figSize=[25,10]','poster=1');

savepdf(h,fullfile('~/proj/afcom/figures/','stimulus.pdf'));