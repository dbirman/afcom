%% BASIC ANALYSIS: AFCOM
% Load the data and do some basic analyses for each subject.
% (1) Check for bias in their response distribution.
% (2) Fit a von Mises to each of the five conditions and compare these


figure(1); clf;
figure(2); clf;
figure(3); clf;

cmap_ = colorblindmap/255;

cmap(5,:) = [0.5 0.5 0.5];
cmap(4,:) = cmap_(7,:);
cmap(3,:) = cmap_(4,:);
cmap(2,:) = cmap_(8,:);
cmap(1,:) = cmap_(1,:);

alldata = [];

for si = 1:length(subjects)
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
        mu = median(dist);
        ci = bootci(100,@median,dist);
        plot(gi-1,dist,'o','MarkerSize',8,'MarkerFaceColor',[0.75 0.75 0.75],'MarkerEdgeColor','w');
        errorbar(gi-1,mu,ci(2)-mu,'-k');
        plot(gi-1,mu,'o','MarkerSize',15,'MarkerFaceColor','k','MarkerEdgeColor','w');
        ticks{end+1} = sprintf('%i-%i',colorgroups(gi-1),colorgroups(gi));
    end
    a = axis;
    axis([a(1) a(2) 0 2]);
    set(gca,'XTick',1:(length(colorgroups)-1),'XTickLabel',ticks);
    
    
    subplot(length(subjects),2,(si-1)*2+2); hold on
    ticks = {};
    for gi = 2:length(dirgroups)
        dat = dirdata(dirgroups(gi-1):dirgroups(gi),:);
        dist = angdist(dat(:,12),dat(:,6));
        mu = median(dist);
        ci = bootci(100,@median,dist);
        plot(gi-1,dist,'o','MarkerSize',8,'MarkerFaceColor',[0.75 0.75 0.75],'MarkerEdgeColor','w');
        errorbar(gi-1,mu,ci(2)-mu,'-k');
        plot(gi-1,mu,'o','MarkerSize',15,'MarkerFaceColor','k','MarkerEdgeColor','w');
        ticks{end+1} = sprintf('%i-%i',dirgroups(gi-1),dirgroups(gi));
    end
    a = axis;
    axis([a(1) a(2) 0 2]);
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
    
    %% Separate color and direction
    colordata = sel(adata,3,1); % report color
    dirdata = sel(adata,3,2);
    
    %% Check for RT differences
    for ti = 1:5
        tdata = colordata(colordata(:,2)==(ti-1),16);
        mu = nanmean(tdata);
        ci = bootci(1000,@nanmean,tdata);
        disp(sprintf('Color RT: %1.2f [%1.2f, %1.2f]',mu,ci(1),ci(2)));
    end
    for ti = 1:5
        tdata = dirdata(dirdata(:,2)==(ti-1),16);
        mu = nanmean(tdata);
        ci = bootci(1000,@nanmean,tdata);
        disp(sprintf('Direction RT: %1.2f [%1.2f, %1.2f]',mu,ci(1),ci(2)));
    end


    %% Fit model
    colordata = sel(adata,3,1); % report color
    dirdata = sel(adata,3,2);
    
    cues = [1 2];
    reportType = {'color','direction'};
    durations = [1 0.25];
    durationType = {'easy','hard'};
    for ci = 1:length(cues)
        for di = 1:length(durations)
            data = sel(adata,3,cues(ci));
            data = sel(data,4,durations(di));
            fit{ci,di} = ac_fitVonMises(data,'bads,nocv');%'nocv');
            fit{ci,di}.dataType = sprintf('%s report %s',durationType{di},reportType{ci});
        end
    end
    
    fits{si} = fit;

    %% Figure for each subject
    figure(2);
    for ci = 1:length(cues)
        for di = 1:length(durations)
            subplot(length(subjects)+1,4,(si-1)*4+(ci-1)*2+di); % plot 1/2 are hard/easy COLOR, then DIRECTION
            hold on
            for tt = 1:5
                if any(tt==[1 4 5])
                    p(tt) = plot(fit{ci,di}.x,fit{ci,di}.out(tt,:),'Color',cmap(tt,:),'LineWidth',1);
                else
                    p(tt) = plot(fit{ci,di}.x,fit{ci,di}.out(tt,:),'Color',cmap(tt,:),'LineWidth',2);
                end
            end
            axis([-pi pi 0 1.5]);
            if ci==1 && di==1
                ylabel('PDF (a.u.)');
            end
            if si==1
                title(sprintf('%s report %s',durationType{di},reportType{ci}));
            end
        end
    end
    legend(fliplr(p),fliplr(fit{1,1}.trialTypes));
end

%% ANALYSIS OF COMBINED DATA

for ci = 1:length(cues)
    for di = 1:length(durations)
        data = sel(alldata,3,cues(ci));
        data = sel(data,4,durations(di));
        fit{ci,di} = ac_fitVonMises(data,'nocv');
        fit{ci,di}.dataType = sprintf('%s report %s',durationType{di},reportType{ci});
    end
end


%% Figure for all subject

figure(2);
for ci = 1:length(cues)
    for di = 1:length(durations)
        subplot(length(subjects)+1,4,length(subjects)*4+(ci-1)*2+di); % plot 1/2 are hard/easy COLOR, then DIRECTION
        hold on
        for tt = 1:5
            if any(tt==[1 4 5])
                p(tt) = plot(fit{ci,di}.x,fit{ci,di}.out(tt,:),'Color',cmap(tt,:),'LineWidth',1);
            else
                p(tt) = plot(fit{ci,di}.x,fit{ci,di}.out(tt,:),'Color',cmap(tt,:),'LineWidth',2);
            end
        end
        axis([-pi pi 0 1.5]);
        xlabel('Distance from target (deg)');
        if ci==1 && di==1
            ylabel('PDF (a.u.)');
        end
    end
end
legend(fliplr(p),fliplr(fit{1,1}.trialTypes));


%%
h = figure(1);
savepdf(h,fullfile('~/proj/afcom/figures/basic_hist.pdf'));
h = figure(2);
savepdf(h,fullfile('~/proj/afcom/figures/basic_VM.pdf'));
h = figure(3);
savepdf(h,fullfile('~/proj/afcom/figures/basic_learn.pdf'));

%% Generate subject data fits
for si = 1:length(subjects)
    h(si) = ac_plotADATA(headers,adatas{si},fits{si});
    savepdf(h(si),fullfile('~/proj/afcom/figures/',sprintf('subj%i_basic.pdf',si)));
end