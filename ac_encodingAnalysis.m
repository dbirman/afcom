%% BASIC ANALYSIS: AFCOM
% Load the data and do some basic analyses for each subject.
% (1) Check for bias in their response distribution.
% (2) Fit a von Mises to each of the five conditions and compare these


figure(1); clf;

cmap_ = colorblindmap/255;

cmap(5,:) = [0.5 0.5 0.5];
cmap(4,:) = cmap_(7,:);
cmap(3,:) = cmap_(4,:);
cmap(2,:) = cmap_(8,:);
cmap(1,:) = cmap_(1,:);

alldata = [];
adatas = {};

infos = {};
for si = 1:length(subjects)
    %% Load data
    [headers,adata,amt] = ac_loadBehavioralData(SIDs{si});
   
    %% Remove dead trials
    adata = adata(~adata(:,5),:);
    
    %% Save
    adatas{si} = adata;
    
    %% Concatenate w/ all data
    if size(adata,1)>400
        alldata = [alldata ; adata];
    end
    
    %% Display count
    disp(sprintf('Subject %s: %i trials',SIDs{si},size(adata,1)));
    
    %% Setup model fits
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
            info = struct;
            info.data = data;
            info.call = 'nocv,bdist';
            info.ci = ci;
            info.di = di;
            info.subj = si;
            info.dataType = sprintf('%s report %s',durationType{di},reportType{ci});
            infos{end+1} = info;
%             fit{ci,di} = ac_fitEncodingModel(data,'nocv,bdist');%'nocv');
            disp(sprintf('%s: %i trials',info.dataType,size(data,1)));
        end
    end
end

%% Add combined data fits
% 
for ci = 1:length(cues)
    for di = 1:length(durations)
        data = sel(alldata,3,cues(ci));
        data = sel(data,4,durations(di));
        info = struct;
        info.data = data;
        info.call = 'nocv,bdist';
        info.ci = ci;
        info.di = di;
        info.subj = -1; % all subjects
        info.dataType = sprintf('%s report %s',durationType{di},reportType{ci});
        infos{end+1} = info;
    end
end

%% Do all the fits
disp(sprintf('Running fits for %i runs',length(infos)));
% sort by the number of trials in each one
len = zeros(size(infos));
for ii = 1:length(infos)
    len(ii) = size(infos{ii}.data,1);
end
[~,idxs] = sort(len,'descend');
infos = infos(idxs);
% run the actual fitting procedure
parfor ii = 1:length(infos)
    infos{ii}.fit = ac_fitEncodingModel(infos{ii}.data,infos{ii}.call);
    infos{ii}.fit.dataType = infos{ii}.dataType;
end
disp('Fits complete');

%% Check that the fits 

%% Sort all the fits

allfits = {}; fits = {};
for si = [-1 1:length(subjects)]
    fit = {};
    for ci = 1:length(cues)
        for di = 1:length(durations)
            for ii = 1:length(infos)
                if infos{ii}.ci==ci && infos{ii}.di==di && infos{ii}.subj==si
                    if si==-1
                        allfit{ci,di} = infos{ii}.fit;
                    else
                        fit{ci,di} = infos{ii}.fit;
                    end
                end
            end
        end
    end
    if si>-1
        fits{si} = fit;
    end
end

%% Plot the fits
% 
% 
for si = 1:length(subjects)
    %% Figure for each subject
    figure(1);
    fit = fits{si};
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
            set(gca,'XTick',[-pi -pi/2 0 pi/2 pi],'XTickLabel',{'-pi','-pi/2','0','pi/2','pi'});
            vline(-pi/2,'--k');
            text(-pi/2,1.25,'side');
            vline(pi*1/3,'--k');
            vline(pi*2/3,'--k');
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
% 
% 
%% Figure for all subject

figure(1);
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
        set(gca,'XTick',[-pi -pi/2 0 pi/2 pi],'XTickLabel',{'-pi','-pi/2','0','pi/2','pi'});
        vline(-pi/2,'--k');
        text(-pi/2,1.25,'side');
        vline(pi*1/3,'--k');
        vline(pi*2/3,'--k');
        xlabel('Distance from target (deg)');
        if ci==1 && di==1
            ylabel('PDF (a.u.)');
        end
    end
end
legend(fliplr(p),fliplr(fit{1,1}.trialTypes));
% 
% 
% %%
h = figure(1);
savepdf(h,fullfile('~/proj/afcom/figures/encoding_VM.pdf'));
% h = figure(2);
% savepdf(h,fullfile('~/proj/afcom/figures/basic_VM.pdf'));
% h = figure(3);
% savepdf(h,fullfile('~/proj/afcom/figures/basic_learn.pdf'));
