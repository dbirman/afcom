%% BASIC ANALYSIS: AFCOM_AVG
% Load the data and do some basic analyses for each subject.
% (1) Check for bias in their response distribution.
% (2) Fit a von Mises to each of the five conditions and compare these

aca_setup;


addpath(genpath('~/proj/afcom'));

cmap_ = colorblindmap/255;

alldata = [];

adatas = {};

for si = 1:length(subjects)
    %% Load data
    mglSetSID(subjects(si));
    [headers,adata] = aca_loadBehavioralData();
   
    %% Remove dead trials
    adata = adata(~adata(:,2),:);
    
    adatas{si} = adata;
    
    l(si) = size(adata,1);
    
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
cis = bootci(10000,@nanmean,acs);
cif = bootci(10000,@nanmean,acf);
%% Figure for all subject
cmap = colorblindmap/255;

poffset = 0.005;
noffset = 0.005;

h = figure(1); clf; hold on

errbar(pscale(xs)+poffset,mean(acs),cis(2,:)-mean(acs),'-','Color',cmap(2,:));
ps(1) = plot(pscale(xs)+poffset,mean(acs),'-','Color',cmap(2,:),'MarkerEdgeColor','w','MarkerSize',5);

errbar(pscale(xs)-noffset,mean(acf),cif(2,:)-mean(acf),'-','Color',cmap(3,:));
ps(2) = plot(pscale(xs)-noffset,mean(acf),'-','Color',cmap(3,:),'MarkerEdgeColor','w','MarkerSize',5);

ylabel('Probability density (a.u.)');
xlabel('Response distance from target (normalized psychological distance)');
set(gca,'XTick',0:0.2:1);
set(gca,'YTick',0:.1:.3);

legend(ps,{'Cue side','Cue color'});
axis([0 1 -0.01 0.3]);

drawPublishAxis('figSize=[8.9,4.5]','poster=0','xLabelOffset=-6/64');

savepdf(h,fullfile('~/proj/afcom/figures','aca.pdf'));

%% plot the circular version of the above data
cmap = colorblindmap/255;

plot_tcc_data(xs,acs,cmap(2,:),'aca_spatial.pdf');
plot_tcc_data(xs,acf,cmap(3,:),'aca_feature.pdf');

%% Test
figure;
hold on
xs = 0:pi/64:pi;
rads = ds(:,4);

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

dbins = linspace(0.20,0.75,3);
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
    plot(pscale(x),n,'-','Color',cmap(di-1,:),'MarkerEdgeColor','w');
%     title(sprintf('%1.2f, dprime %1.2f',mid,dp));
end
a = axis;
axis([0 1 0 a(4)]);
legend({'0.25 - 0.5 s','0.5 - 0.75 s'});

ylabel('Probability density (a.u.)');
xlabel('Psychophysical distance from target (normalized)');
set(gca,'XTick',0:1,'XTickLabel',0:1);
set(gca,'YTick',0:.1:.2);

drawPublishAxis('figSize=[8.9,4.5]','poster=0');

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
    plot(pscale(x),n,'-','Color',cmap(di-1,:),'MarkerEdgeColor','w');
%     title(sprintf('%1.2f, dprime %1.2f',mid,dp));
end
a = axis;
axis([0 1 0 a(4)]);
legend({'0.25 - 0.5 s','0.5 - 0.75 s'});

ylabel('Probability density (a.u.)');
xlabel('Psychophysical distance from target (normalized)');
set(gca,'XTick',0:1,'XTickLabel',[0 1]);
set(gca,'YTick',0:.1:.2);
drawPublishAxis('figSize=[8.9,4.5]','poster=0');

savepdf(h,fullfile('~/proj/afcom/figures','distance.pdf'));

%%
%%%%%%% FINISHED
return
%%

% I don't fit the full model anymore to the ACA data. Just the raw fits.
% Note that there's no "no attention" condition here, so the best we can do
% is just compare these directly to each other. 

% 
% %% Test the full model fit
% fit = aca_fitTCCModel(adata,'nocv,bads',[]);
% aca_plotTCCModel(fit);
% 
% 
% %% Test the population variant model fit
% fit = aca_fitTCCPopulationModel(adata,'nocv,bads',[]);
% % aca_plotTCCModel(fit);
% 
% %% test output from fit
% figure(1); clf; hold on
% xs = -pi:pi/128:pi;
% ys = computeTCCfromPopulation(xs,fitpop.params.sigma,fitpop.params.p);
% plot(xs,ys./sum(ys),'-b')
% ys2 = computeTCCPDF(xs,fit.params.p);
% plot(xs,ys2,'-r')