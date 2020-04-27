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
clear dprimes dprimef acs acf dprimes_target dprimef_target
for ai = 1:length(adatas)
    ds = adatas{ai}(adatas{ai}(:,3)==1,:);
    df = adatas{ai}(adatas{ai}(:,3)==2,:);
    
    % fit the TCC model
    dprimes(ai) = fitTCC(ds(:,4));
    dprimef(ai) = fitTCC(df(:,4));
    
    % also fit TCC separately for blue/yellow and left/right
    for target = 1:2
        ds_t = ds(ds(:,8)==target,:);
        dprimes_target(ai,target) = fitTCC(ds_t(:,4));
        df_t = df(df(:,8)==target,:);
        dprimef_target(ai,target) = fitTCC(df_t(:,4));
    end
    
%     kappas(ai) = fitVonMisesEstimation(ds(:,4));
%     kappaf(ai) = fitVonMisesEstimation(df(:,4));
    
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

%% Make a figure showing how well the fits worked (check that the model is actuall working?)
h = figure; hold on

plot(xs,mean(acs),'o','MarkerFaceColor','k','MarkerEdgeColor','w');
plot(xs,computeTCCPDF(xs,mean(dprimes)));
y = vonMises(xs,0,mean(kappas));
y = y./sum(y);
plot(xs,y);

%% info for text
cidprimes = bootci(10000,@nanmean,dprimes);
cidprimef = bootci(10000,@nanmean,dprimef);
disp(sprintf('Mean spatial d'' value: %1.2f 95%% CI [%1.2f, %1.2f]',mean(dprimes),cidprimes(1),cidprimes(2)));
disp(sprintf('Mean feature d'' value: %1.2f 95%% CI [%1.2f, %1.2f]',mean(dprimef),cidprimef(1),cidprimef(2)));

ci = bootci(10000,@nanmean,dprimes-dprimef);

disp(sprintf('Difference in d'': %1.2f 95%% CI [%1.2f, %1.2f]',mean(dprimes-dprimef),ci(1), ci(2)));
%% Figure for all subject
cmap = colorblindmap/255;

poffset = 0.005;
noffset = 0.005;


h = figure(1); clf; hold on

ex = pscale(xs)+poffset;
ex = [ex fliplr(ex)];
h2 = fill(ex,[cis(1,:) fliplr(cis(2,:))],cmap(2,:));
h2.FaceAlpha = 0.25;
h2.LineStyle = 'none';
% errbar(pscale(xs)+poffset,mean(acs),cis(2,:)-mean(acs),'-','Color',cmap(2,:));
ps(1) = plot(pscale(xs)+poffset,mean(acs),'-','Color',cmap(2,:),'MarkerEdgeColor','w','MarkerSize',5);
% y = computeTCCPDF(xs,mean(dprimes));
% plot(pscale(xs*180/pi),y,'-','Color',cmap(2,:));

ex = pscale(xs)+poffset;
ex = [ex fliplr(ex)];
h2 = fill(ex,[cif(1,:) fliplr(cif(2,:))],cmap(3,:));
h2.FaceAlpha = 0.25;
h2.LineStyle = 'none';
% errbar(pscale(xs)-noffset,mean(acf),cif(2,:)-mean(acf),'-','Color',cmap(3,:));
ps(2) = plot(pscale(xs)-noffset,mean(acf),'-','Color',cmap(3,:),'MarkerEdgeColor','w','MarkerSize',5);
% y = computeTCCPDF(xs,mean(dprimef));
% plot(pscale(xs*180/pi),y,'-','Color',cmap(3,:));

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






%% DURATION PLOTS

%% First separate the data by subject and get each subject's histogram

xs = pi/64:pi/32:pi;
dbins = linspace(0.20,0.75,3);

dur_out = nan(length(adatas),2,length(xs));
dp = zeros(length(adatas),2);
for ai = 1:length(adatas)
    data = adatas{ai};
    
    for di = 2:length(dbins)
        low = dbins(di-1);
        high = dbins(di);
        mid = mean([low high]);
        idxs = (data(:,5)>=low).*(data(:,5)<high);
        dat = data(logical(idxs),:);
        [n,x] = hist(dat(:,4),xs);
        n = n./sum(n);
        
        dur_out(ai,di-1,:) = n;
        
        dp(ai,di-1) = fitTCC(dat(:,4));
    end
end
% now average the histograms
dur_ = squeeze(mean(dur_out));
dur_ci = squeeze(bootci(1000,@nanmean,dur_out));

% plot
h = figure; hold on

cmap = brewermap(13,'Purples');
cmap = cmap([7 13],:);



for di = 1:2
    % error bars
    ex = pscale(x);
    ex = [ex fliplr(ex)];
    cis = [squeeze(dur_ci(1,di,:))' fliplr(squeeze(dur_ci(2,di,:))')];
    h2 = fill(ex,cis,cmap(di,:));
    h2.FaceAlpha = 0.25;
    h2.LineStyle = 'none';
    % lines
    p(di) = plot(pscale(x),dur_(di,:),'-','Color',cmap(di,:),'MarkerEdgeColor','w');
end


a = axis;
axis([0 1 0 a(4)]);
legend(p,{'0.25 - 0.5 s','0.5 - 0.75 s'});

ylabel('Probability density (a.u.)');
xlabel('Psychophysical distance from target (normalized)');
set(gca,'XTick',0:1,'XTickLabel',0:1);
set(gca,'YTick',0:.1:.2);

drawPublishAxis('figSize=[4.5,4.5]','poster=0');

% % savepdf(h,fullfile('~/proj/afcom/figures','duration.pdf'));


%% info for text

ddur = dp(:,2)-dp(:,1);
cidur = bootci(10000,@nanmean,ddur);
disp(sprintf('Difference in d'' due to duration: %1.2f 95%% CI [%1.2f, %1.2f]',mean(ddur),cidur(1),cidur(2)));

% signed-rank test
p = signrank(dp(:,1),dp(:,2));

disp(sprintf('signed-rank %P=%1.2f',p));

%% DISTANCE PLOTS

%% First separate the data by subject and get each subject's histogram

xs = pi/64:pi/32:pi;
dbins = linspace(0,0.75*pi,3);

dp = zeros(length(adatas),2);
dur_out = nan(length(adatas),2,length(xs));
for ai = 1:length(adatas)
    data = adatas{ai};
    dist = angdist(data(:,6),data(:,7));
    
    for di = 2:length(dbins)
        low = dbins(di-1);
        high = dbins(di);
        mid = mean([low high]);
        idxs = (dist>=low).*(dist<high);
        dat = data(logical(idxs),:);
        [n,x] = hist(dat(:,4),xs);
        n = n./sum(n);
        
        dur_out(ai,di-1,:) = n;
        
        dp(ai,di-1) = fitTCC(dat(:,4));
    end
end
% now average the histograms
dur_ = squeeze(mean(dur_out));
dur_ci = squeeze(bootci(1000,@nanmean,dur_out));

%% plot
h = figure; hold on

cmap = brewermap(13,'Oranges');
cmap = cmap([7 13],:);

for di = 1:2
    % error bars
    ex = pscale(x);
    ex = [ex fliplr(ex)];
    cis = [squeeze(dur_ci(1,di,:))' fliplr(squeeze(dur_ci(2,di,:))')];
    h2 = fill(ex,cis,cmap(di,:));
    h2.FaceAlpha = 0.25;
    h2.LineStyle = 'none';
    % lines
    p(di) = plot(pscale(x),dur_(di,:),'-','Color',cmap(di,:),'MarkerEdgeColor','w');
end


a = axis;
axis([0 1 0 0.25]);
legend(p,{'0 - pi/3','pi/3 - pi'});

ylabel('Probability density (a.u.)');
xlabel('Psychophysical distance from target (normalized)');
set(gca,'XTick',0:1,'XTickLabel',0:1);
set(gca,'YTick',0:.1:.2);

drawPublishAxis('figSize=[4.5,4.5]','poster=0');

savepdf(h,fullfile('~/proj/afcom/figures','distance.pdf'));

%%
ddur = dp(:,1)-dp(:,2);
cidur = bootci(10000,@nanmean,ddur);
disp(sprintf('Difference in d'' due to distance: %1.2f 95%% CI [%1.2f, %1.2f]',mean(ddur),cidur(1),cidur(2)));

p = signrank(dp(:,1),dp(:,2));

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