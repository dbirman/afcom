%% BASIC ANALYSIS: AFCOM_AVG
% Load the data and do some basic analyses for each subject.
% (1) Check for bias in their response distribution.
% (2) Fit a von Mises to each of the five conditions and compare these
addpath(genpath('~/proj/afcom'));

cmap_ = colorblindmap/255;

cmap(1,:) = cmap_(7,:);
cmap(2,:) = cmap_(4,:);

alldata = [];

adatas = {};

for si = 1:length(subjects)
    %% Load data
    mglSetSID(subjects(si));
    [headers,adata] = aca_loadBehavioralData();
   
    %% Remove dead trials
    adata = adata(~adata(:,2),:);
    
    adatas{si} = adata;
    
    %% Remove the first two runs of each type
%     runs = unique(adata(:,15));
%     ccount = 0;
%     dcount = 0;
%     for ri = 1:length(runs)
%         if 
%     end
    
    %% Concatenate w/ all data
    if size(adata,1)>100
        alldata = [alldata ; adata];
    end
end

%% Compute mean and SD for the histograms
xs = pi/64:pi/32:pi;
for ai = 1:length(adatas)
    ds = adatas{ai}(adatas{ai}(:,3)==1,:);
    df = adatas{ai}(adatas{ai}(:,3)==2,:);
    
    cs = hist(ds(:,4),xs);
    cs = cs ./ sum(cs);
    
    cf = hist(df(:,4),xs);
    cf = cf ./ sum(cf);
    
    acs(ai,:) = cs;
    acf(ai,:) = cf;
end

% compute errbars
cis = bootci(1000,@nanmean,acs);
cif = bootci(1000,@nanmean,acf);
%% Figure for all subject
cmap = colorblindmap/255;

offset = pi/128;

h = figure; hold on
% subplot(211); hold on
errbar(xs+offset,mean(acs),cis(2,:)-mean(acs),'-','Color',cmap(2,:));
ps(1) = plot(xs+offset,mean(acs),'o','MarkerFaceColor',cmap(2,:),'MarkerEdgeColor','w','MarkerSize',8);

% add the model fit of the TCC model
% dp = fitTCC(dat_spatial(:,4));
% plot(xs,computeTCCPDF(xs,dp),'-k');
% legend({num2str(dp)});
ylabel('Probability density (a.u.)');
xlabel('Response distance from target (degs)');

% vline(pi/2,'--r');


% subplot(212); hold on
errbar(xs,mean(acf),cif(2,:)-mean(acf),'-','Color',cmap(3,:));
ps(2) = plot(xs,mean(acf),'o','MarkerFaceColor',cmap(3,:),'MarkerEdgeColor','w','MarkerSize',8);

% dp = fitTCC(dat_feature(:,4));
% plot(xs,computeTCCPDF(xs,dp),'-k');
% legend({num2str(dp)});
ylabel('Probability density (a.u.)');
xlabel('Response distance from target (degs)');

% v = vline(median(dat_spatial(:,4)),'--');
% set(v,'Color',cmap(2,:));
% v = vline(median(dat_feature(:,4)),'--');
% set(v,'Color',cmap(3,:));
set(gca,'XTick',0:pi/4:pi,'XTickLabel',180/pi*(0:pi/4:pi));
set(gca,'YTick',0:.1:.3);
% vline(pi/2,'--r');
legend(ps,{'Cue side','Cue color'});
axis([0 pi -0.01 0.3]);

drawPublishAxis('figSize=[30,20]','labelFontSize=18','xLabelOffset=-6/64');

savepdf(h,fullfile('~/proj/afcom/figures','aca.pdf'));

%% Test
figure;
hold on
xs = 0:pi/64:pi;
rads = dat_feature(:,4);

[n,x] = hist(rads,xs);
n = n/sum(n);
plot(x,n,'o','MarkerFaceColor','k');
l = [];
for dprime=0.1:.25:2
    plot(xs,computeTCCPDF(xs,dprime));
    l(end+1) = -sum(log(computeTCCPDF(rads,dprime)));
end
legend(l);

%% Block the data by duration quantile and plot 
figure;
xs = 0:pi/64:pi;

dbins = linspace(0.25,0.75,8);
for di = 2:length(dbins)
    subplot(length(dbins)-1,1,di-1); hold on
    
    low = dbins(di-1);
    high = dbins(di);
    mid = mean([low high]);
    idxs = (alldata(:,5)>=low).*(alldata(:,5)<high);
    dat = alldata(logical(idxs),:);
    
    dp = fitTCC(dat(:,4));
    plot(xs,computeTCCPDF(xs,dp),'-k');
    
    [n,x] = hist(dat(:,4),xs);
    n = n./sum(n);
    plot(x,n,'ok');
    axis([0 pi 0 0.25]);
    title(sprintf('%1.2f, dprime %1.2f',mid,dp));
end

%% Block the data by distance quantile and plot
figure;
xs = 0:pi/64:pi;

dist = angdist(alldata(:,6),alldata(:,7));

dbins = linspace(0,0.75*pi,8);
for di = 2:length(dbins)
    subplot(length(dbins)-1,1,di-1); hold on
    
    low = dbins(di-1);
    high = dbins(di);
    mid = mean([low high]);
    idxs = (dist>=low).*(dist<high);
    dat = alldata(logical(idxs),:);
    
    dp = fitTCC(dat(:,4));
    plot(xs,computeTCCPDF(xs,dp),'-k');
    
    [n,x] = hist(dat(:,4),xs);
    n = n./sum(n);
    plot(x,n,'ok');
    axis([0 pi 0 0.25]);
    title(sprintf('%1.2f, dprime %1.2f',mid,dp));
end

%% Test the full model fit
fit = aca_fitTCCModel(adata,'nocv,bads',[]);
aca_plotTCCModel(fit);



%% example plot for poster
cmap = brewermap(15,'RdYlGn');
cmap = shift(cmap,[7,0]);
cmap = cmap(1:13,:);

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
    plot(mu,activation,'o','MarkerFaceColor',cmap(ui,:),'MarkerEdgeColor','w','MarkerSize',5);
end

ylabel('Tuning profile (a.u.)');
xlabel('Distance from presented stimulus (deg)');
axis([-pi pi -0.1 1.1]);
set(gca,'XTick',[-pi 0 pi]);
set(gca,'YTick',[0 1]);
drawPublishAxis('figSize=[50,10]','labelFontSize=18');

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