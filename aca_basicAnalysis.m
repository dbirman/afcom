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

h = figure(1); hold on
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

% v = vline(median(dat_spatial(:,4)),'--');
% set(v,'Color',cmap(2,:));
% v = vline(median(dat_feature(:,4)),'--');
% set(v,'Color',cmap(3,:));
ylabel('Probability density (a.u.)');
xlabel('Response distance from target (degs)');
set(gca,'XTick',0:pi/4:pi,'XTickLabel',180/pi*(0:pi/4:pi));
set(gca,'YTick',0:.1:.3);
% vline(pi/2,'--r');
legend(ps,{'Cue side','Cue color'});
axis([0 pi -0.01 0.3]);

drawPublishAxis('figSize=[40,12]','poster=1');

savepdf(h,fullfile('~/proj/afcom/figures','aca.pdf'));

%% plot the circular version of the above data
cmap = colorblindmap/255;

plot_tcc_data(xs,acs,cmap(2,:),'aca_spatial.pdf');
plot_tcc_data(xs,acf,cmap(3,:),'aca_feature.pdf');

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
h = figure; hold on
xs = pi/64:pi/32:pi;

cmap = brewermap(13,'Purples');
cmap = cmap([7 13],:);

dbins = linspace(0.25,0.75,3);
for di = 2:length(dbins)
    
    low = dbins(di-1);
    high = dbins(di);
    mid = mean([low high]);
    idxs = (alldata(:,5)>=low).*(alldata(:,5)<high);
    dat = alldata(logical(idxs),:);
    
%     dp = fitTCC(dat(:,4));
%     plot(xs,computeTCCPDF(xs,dp),'-k');
    
    [n,x] = hist(dat(:,4),xs);
    n = n./sum(n);
    plot(x,n,'o','MarkerFaceColor',cmap(di-1,:),'MarkerEdgeColor','w');
%     title(sprintf('%1.2f, dprime %1.2f',mid,dp));
end
a = axis;
axis([0 pi 0 a(4)]);
legend({'0.25 - 0.5 s','0.5 - 0.75 s'});

ylabel('Probability density (a.u.)');
xlabel('Response distance from target (degs)');
set(gca,'XTick',0:pi/4:pi,'XTickLabel',180/pi*(0:pi/4:pi));
set(gca,'YTick',0:.1);

drawPublishAxis('figSize=[20,10]','poster=1');

savepdf(h,fullfile('~/proj/afcom/figures','duration.pdf'));

%% Block the data by distance quantile and plot
h = figure; hold on
xs = pi/64:pi/32:pi;

cmap = brewermap(13,'Oranges');
cmap = cmap([7 13],:);
dist = angdist(alldata(:,6),alldata(:,7));

dbins = linspace(0,0.75*pi,3);
for di = 2:length(dbins)
    
    low = dbins(di-1);
    high = dbins(di);
    mid = mean([low high]);
    idxs = (dist>=low).*(dist<high);
    dat = alldata(logical(idxs),:);
    
%     dp = fitTCC(dat(:,4));
%     plot(xs,computeTCCPDF(xs,dp),'-k');
    
    [n,x] = hist(dat(:,4),xs);
    n = n./sum(n);
    plot(x,n,'o','MarkerFaceColor',cmap(di-1,:),'MarkerEdgeColor','w');
%     title(sprintf('%1.2f, dprime %1.2f',mid,dp));
end
a = axis;
axis([0 pi 0 a(4)]);
legend({'0 - 67.5 deg','67.5 - 135 deg'});

ylabel('Probability density (a.u.)');
xlabel('Response distance from target (degs)');
set(gca,'XTick',0:pi/4:pi,'XTickLabel',180/pi*(0:pi/4:pi));
set(gca,'YTick',0:.1);

drawPublishAxis('figSize=[20,10]','poster=1');

savepdf(h,fullfile('~/proj/afcom/figures','distance.pdf'));

%% Test the full model fit
fit = aca_fitTCCModel(adata,'nocv,bads',[]);
aca_plotTCCModel(fit);


%% Test the population variant model fit
fit = aca_fitTCCPopulationModel(adata,'nocv,bads',[]);
% aca_plotTCCModel(fit);