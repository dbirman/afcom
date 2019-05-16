

%% example plot for poster
cmap = brewermap(15,'RdYlGn');
cmap = shift(cmap,[7,0]);
cmap = cmap(1:12,:);

% plot a TCC channel model, showing the activation of different units 
units = -pi:(2*pi/(size(cmap,1)-1)):pi;
units = sort(units);
x = -pi:pi/128:pi;

h = figure(2); clf; hold on
vline(0,'--k');
for ui = 1:length(units)
    mu = units(ui);
    profile = 1-pscale(abs(x-mu)*180/pi);
    % compute my activation at 0
    activation = 1-pscale(abs(mu)*180/pi);
    plot(x,profile,'-','Color',cmap(ui,:));
end
for ui = 1:length(units)
    mu = units(ui);
    profile = 1-pscale(abs(x-mu)*180/pi);
    % compute my activation at 0
    activation = 1-pscale(abs(mu)*180/pi);
    plot([mu 0],[activation activation],'--','Color',cmap(ui,:));
end
for ui = 1:length(units)
    mu = units(ui);
    profile = 1-pscale(abs(x-mu)*180/pi);
    % compute my activation at 0
    activation = 1-pscale(abs(mu)*180/pi);
    plot(mu,activation,'o','MarkerFaceColor',cmap(ui,:),'MarkerEdgeColor','w','MarkerSize',10);
end

ylabel('Tuning profile (a.u.)');
xlabel('Distance from presented stimulus (deg)');
axis([-pi pi -0.1 1.1]);
set(gca,'XTick',[-pi 0 pi]);
set(gca,'YTick',[0 1]);
drawPublishAxis('figSize=[25,10]','poster=1');

savepdf(h,fullfile('~/proj/afcom/figures/','TCC model.pdf'));

%% same example plot, but with just the activations + noise

cmap = brewermap(15,'RdYlGn');

% plot a TCC channel model, showing the activation of different units 
units = -pi:(2*pi/14):pi;
units = sort(units);
x = -pi:pi/128:pi;

h = figure(2); clf; hold on
for ui = 1:length(units)
    mu = units(ui);
    profile = 1-pscale(abs(x-mu)*180/pi);
    % compute my activation at 0
    activation = 1-pscale(abs(mu)*180/pi);
    plot(mu,activation,'o','MarkerFaceColor','k','MarkerEdgeColor','w','MarkerSize',5);
    errbar(mu,activation,0.2,'-k');
end

ylabel('Tuning profile (a.u.)');
xlabel('Distance from presented stimulus (deg)');
axis([-pi pi -0.3 1.3]);
set(gca,'XTick',[-pi 0 pi]);
set(gca,'YTick',[0 1]);
drawPublishAxis('figSize=[5,3]');

savepdf(h,fullfile('~/proj/afcom/figures/','TCC activations.pdf'));


cmap = brewermap(15,'RdYlGn');

%% plot a TCC channel model, showing the activation of different units 
units = -pi:(2*pi/14):pi;
units = sort(units);
x = -pi:pi/128:pi;

h = figure(2); clf; hold on
for ui = 1:length(units)
    mu = units(ui);
    % compute my activation at 0
    activation = 1-pscale(abs(mu)*180/pi);
    plot(mu,activation,'o','MarkerFaceColor','k','MarkerEdgeColor','w','MarkerSize',10);
    errbar(mu,activation,0.2,'-k');
    
    activation = 1-pscale(abs(angdist(mu,0.75*pi))*180/pi);
    plot(mu,activation,'o','MarkerFaceColor','r','MarkerEdgeColor','w','MarkerSize',10);
    errbar(mu,activation,0.2,'-r');
end

ylabel('Tuning profile (a.u.)');
xlabel('Distance from presented stimulus (deg)');
axis([-pi pi -0.3 1.3]);
set(gca,'XTick',[-pi 0 pi]);
set(gca,'YTick',[0 1]);
drawPublishAxis('figSize=[25,15]','poster=1');

savepdf(h,fullfile('~/proj/afcom/figures/','TCC activations.pdf'));


%% TCC Channel model with variable sensitivity 
units = -pi:(2*pi/14):pi;
units = sort(units);
x = -pi:pi/128:pi;

h = figure(2); clf; hold on
for ui = 1:length(units)
    mu = units(ui);
    % compute my activation at 0
    activation = 1-pscale(abs(mu)*180/pi);
    errbar(mu,activation,0.2,'-k');
    plot(mu,activation,'o','MarkerFaceColor','k','MarkerEdgeColor','w','MarkerSize',10);
    
    activation = 1-pscale(abs(mu)*180/pi);
    errbar(mu+0.03,1.5*activation,0.2,'-r');
    plot(mu+0.03,1.5*activation,'o','MarkerFaceColor','r','MarkerEdgeColor','w','MarkerSize',10);
end

ylabel('Tuning profile (a.u.)');
xlabel('Distance from presented stimulus (deg)');
axis([-pi pi -0.3 2]);
set(gca,'XTick',[-pi 0 pi],'XTickLabel',[-180 0 180]);
set(gca,'YTick',[0 1]);
drawPublishAxis('figSize=[25,15]','poster=1');

savepdf(h,fullfile('~/proj/afcom/figures/','TCC activations_sensitivity.pdf'));


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
    ax1 = ax-d*cos(angs(ai));
    ax2 = ax+d*cos(angs(ai));
    ay1 = ay-d*sin(angs(ai));
    ay2 = ay+d*sin(angs(ai));
    arrow([ax1 ay1],[ax2,ay2],'Color',cmap.(types{ai}));
end
for ai = 1:length(angs)
    activation(ai,:) = 1-pscale(abs(units-angs(ai))*180/pi);
end
plot(x,activation' * weights','-k');
ylabel('Channel activation (a.u.)');
xlabel('Distance from presented stimulus (deg)');
axis([-pi pi -0.3 1]);
set(gca,'XTick',[-pi 0 pi],'XTickLabel',[-180 0 180]);
set(gca,'YTick',[0 1]);
drawPublishAxis('figSize=[25,15]','poster=1');

savepdf(h,fullfile('~/proj/afcom/figures/','afcom_tcc_activationmodel.pdf'));