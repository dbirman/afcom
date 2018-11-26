%% BASIC ANALYSIS: AFCOM
% Load the data and do some basic analyses for each subject.
% (1) Check for bias in their response distribution.
% (2) Fit a von Mises to each of the five conditions and compare these


figure(1); clf;
figure(2); clf;
figure(3); clf;

cmap = brewermap(5,'Dark2');

alldata = [];

for si = 1%:length(subjects)
    %% Load data
    [headers,adata,amt] = ac_loadBehavioralData(SIDs{si});
   
    %% Remove dead trials
    adata = adata(~adata(:,5),:);
    
    %% Remove the first two runs of each type
%     runs = unique(adata(:,15));
%     ccount = 0;
%     dcount = 0;
%     for ri = 1:length(runs)
%         if 
%     end
    
    %% Concatenate w/ all data
    if size(adata,1)>400
        alldata = [alldata ; adata];
    end
    
    %% Display count
    disp(sprintf('Subject %s: %i trials',SIDs{si},size(adata,1)));
    
    %% Improvement figure
    % split data into color and direction
    colordata = sel(adata,3,1);
    dirdata = sel(adata,3,2);
    % take only target trials trialType==3
    colordata = sel(colordata,2,3);
    dirdata = sel(dirdata,2,3);
    
    % for every 20 trials, plot the abs distribution of respAngle - targetAngle
    colorgroups = 1:20:size(colordata,1);
    dirgroups = 1:20:size(dirdata,1);
    figure(3);
    subplot(length(subjects),2,(si-1)*2+1); hold on
    ticks = {};
    for gi = 2:length(colorgroups)
        dat = colordata(colorgroups(gi-1):colorgroups(gi),:);
        dist = angdist(dat(:,12),dat(:,6));
        mu = mean(dist);
        ci = bootci(100,@mean,dist);
        plot(gi-1,dist,'o','MarkerSize',8,'MarkerFaceColor',[0.75 0.75 0.75],'MarkerEdgeColor','w');
        errorbar(gi-1,mu,ci(2)-mu,'-k');
        plot(gi-1,mu,'o','MarkerSize',15,'MarkerFaceColor','k','MarkerEdgeColor','w');
        ticks{end+1} = sprintf('%i-%i',colorgroups(gi-1),colorgroups(gi));
    end
    set(gca,'XTick',1:(length(colorgroups)-1),'XTickLabel',ticks);
    
    
    subplot(length(subjects),2,(si-1)*2+2); hold on
    ticks = {};
    for gi = 2:length(dirgroups)
        dat = dirdata(dirgroups(gi-1):dirgroups(gi),:);
        dist = angdist(dat(:,12),dat(:,6));
        mu = mean(dist);
        ci = bootci(100,@mean,dist);
        plot(gi-1,dist,'o','MarkerSize',8,'MarkerFaceColor',[0.75 0.75 0.75],'MarkerEdgeColor','w');
        errorbar(gi-1,mu,ci(2)-mu,'-k');
        plot(gi-1,mu,'o','MarkerSize',15,'MarkerFaceColor','k','MarkerEdgeColor','w');
        ticks{end+1} = sprintf('%i-%i',dirgroups(gi-1),dirgroups(gi));
    end
    set(gca,'XTick',1:(length(dirgroups)-1),'XTickLabel',ticks);
    
    %% Distribution figure
    colordata = sel(adata,3,1); % report color
    dirdata = sel(adata,3,2);
    figure(1);
    subplot(length(subjects)+1,2,(si-1)*2+1);
    hist(mod(colordata(:,12),2*pi),16,'FaceColor',[0.75 0.75 0.75],'EdgeColor','w');
    a = axis;
    axis([0 2*pi a(3) a(4)]);
    set(gca,'XTick',0:pi/4:2*pi);
    subplot(length(subjects)+1,2,(si-1)*2+2);
    hist(mod(dirdata(:,12),2*pi),16,'FaceColor',[0.75 0.75 0.75],'EdgeColor','w');
    a = axis;
    axis([0 2*pi a(3) a(4)]);
    set(gca,'XTick',0:pi/4:2*pi);

    %% Fit model
    colordata = sel(adata,3,1); % report color
    dirdata = sel(adata,3,2);
    
    colorfit = ac_fitVonMises(colordata,'lapseall');
    dirfit = ac_fitVonMises(dirdata,'lapseall');
%     %% Group trials by type
%     tType = {'all','spatial','feature','target','baseline'};
%     tData = cell(1,5);
%     for ti = 0:4
%         tData{ti+1} = adata(adata(:,2)==ti,:);
%     end
%     
%     %% Check for RT differences
%     for ti = 1:5
%         mu = nanmean(tData{ti}(:,16));
%         ci = bootci(1000,@nanmean,tData{ti}(:,16));
%         disp(sprintf('RT: %1.2f [%1.2f, %1.2f]',mu,ci(1),ci(2)));
%     end
%     
%     %% Fit a vonMises
%     clear cfit dfit
%     for ti = 1:5
%         disp(sprintf('Fitting VM for trials: %s',tType{ti}));
%         cdat = tData{ti};
%         % get color task
%         colordata = cdat(cdat(:,3)==1,:);
%         dirdata = cdat(cdat(:,3)==2,:);
%         dur = [1 0.25];
%         for di = 1:length(dur)
%             cfit{ti,di} = ac_fitVonMises(colordata(colordata(:,4)==dur(di),:));
%             dfit{ti,di} = ac_fitVonMises(dirdata(dirdata(:,4)==dur(di),:));
%         end
%     end
%     
%     %% Print out the SD estimates
%     disp('Report: direction');
%     for ti = 1:5
%         disp(sprintf('%s SD: %1.2f',tType{ti},sqrt(1/cfit{ti}.params.kappa)));
%     end
%     disp('Report: color');
%     for ti = 1:5
%         disp(sprintf('%s SD: %1.2f',tType{ti},sqrt(1/dfit{ti}.params.kappa)));
%     end
% 
%     %% Figure for each subject
%     figure(2);
%     for di = 1:length(dur)
%         subplot(length(subjects)+1,4,(si-1)*4+0+di);
%         hold on
%         for ti = 1:5
%             plot(cfit{ti,di}.x,cfit{ti,di}.out,'Color',cmap(ti,:));
%         end
%         axis([-pi pi 0 2.5]);
%     %     drawPublishAxis;
%         subplot(length(subjects)+1,4,(si-1)*4+2+di);
%         hold on
%         for ti = 1:5
%             plot(dfit{ti,di}.x,dfit{ti,di}.out,'Color',cmap(ti,:));
%         end
%         a = axis;
%         axis([-pi pi 0 2.5]);
%         legend(tType);
%     %     drawPublishAxis;
%     end
    
end

%% Group trials by type
tType = {'all','spatial','feature','target','baseline'};
tData = cell(1,5);
for ti = 0:4
    tData{ti+1} = alldata(alldata(:,2)==ti,:);
end

%% Fit a vonMises
% clear fit
% for ti = 1:5
%     disp(sprintf('Fitting VM for trials: %s',tType{ti}));
%     cdat = tData{ti};
%     cfit{ti} = ac_fitVonMises(cdat(cdat(:,3)==1,:));
%     dfit{ti} = ac_fitVonMises(cdat(cdat(:,3)==2,:));
% end

clear cfit dfit
for ti = 1:5
    disp(sprintf('Fitting VM for trials: %s',tType{ti}));
    cdat = tData{ti};
    % get color task
    colordata = cdat(cdat(:,3)==1,:);
    dirdata = cdat(cdat(:,3)==2,:);
    dur = [1 0.25];
    for di = 1:length(dur)
        cfit{ti,di} = ac_fitVonMises(colordata(colordata(:,4)==dur(di),:));
        dfit{ti,di} = ac_fitVonMises(dirdata(dirdata(:,4)==dur(di),:));
    end
end

%% Figure for all subject

figure(2);
for di = 1:length(dur)
    subplot(length(subjects)+1,4,length(subjects)*4+0+di);
    hold on
    for ti = 1:5
        plot(cfit{ti,di}.x,cfit{ti,di}.out,'Color',cmap(ti,:));
    end
    axis([-pi pi 0 2.5]);
%     drawPublishAxis;
    subplot(length(subjects)+1,4,length(subjects)*4+2+di);
    hold on
    for ti = 1:5
        plot(dfit{ti,di}.x,dfit{ti,di}.out,'Color',cmap(ti,:));
    end
    a = axis;
    axis([-pi pi 0 2.5]);
    legend(tType);
%     drawPublishAxis;
end


%%
h = figure(2);
savepdf(h,fullfile('~/proj/afcom/figures/basic.pdf'));