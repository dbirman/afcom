%% BASIC ANALYSIS: AFCOM
% Load the data and do some basic analyses for each subject.
% (1) Check for bias in their response distribution.
% (2) Fit a von Mises to each of the five conditions and compare these
ac_setup;

cues = [1 2];

alldata = [];
adatas = {};
infos = {};

for si = 1:length(subjects)
    %% Load data
    [headers,adata] = ac_loadBehavioralData(SIDs{si});
    
% avars = {'target','trialType','cue','duration','dead','targetAngle','distractorAngle','angle1','angle2','angle3',...
%     'angle4','respAngle','respDistance','distDistance'};

% TARGET TRIALTYPE CUE DURATION DEAD TARGETANGLE DISTRACTORANGLE
%    1      2       3     4       5       6           7
% ANGLE1 ANGLE2 ANGLE3 ANGLE4 RESPANGLE RESPDISTANCE DISTDISTANCE RUNS RT
%    8      9     10     11      12          13            14      15  16
   
    %% Remove dead trials
    disp(sprintf('Dropping %i trials that included eye movements',sum(adata(:,5))));
    adata = adata(~adata(:,5),:);
    
    %% Remove any trials with duration > 0.3 (these were "easy" trials in the initial versions)
    disp(sprintf('Dropping %i trials from training or easy mode',sum(adata(:,4)>0.3)));
    adata = adata(adata(:,4)<=0.3,:);
    
    %% Save
    adatas{si} = adata;
    l(si) = size(adata,1);
    
    %% Concatenate w/ all data
    if size(adata,1)>250
        disp(sprintf('Adding data to alldata'));
        alldata = [alldata ; adata];
    end
    
    %% Display count
    disp(sprintf('Subject %s: %i trials',SIDs{si},size(adata,1)));
    
    %% Setup model fits
    colordata = sel(adata,3,1); % report color
    dirdata = sel(adata,3,2);
    
    reportType = {'color','direction'};
    durations = [0.25 0.3];
    durationType = {'hard1','hard2'};
    for ci = 1:length(cues)
        for di = 1:length(durations)
            data = sel(adata,3,cues(ci));
            data = sel(data,4,durations(di));
            if size(data,1)>0
                info = struct;
                info.data = data;
                info.call = 'nocv,bads';
                info.ci = ci;
                info.di = di;
                info.subj = si;
                info.dataType = sprintf('%s report %s',durationType{di},reportType{ci});
                infos{end+1} = info;
    %             fit{ci,di} = ac_fitEncodingModel(data,'nocv,bdist');%'nocv');
                disp(sprintf('%s: %i trials',info.dataType,size(data,1)));
            else
                disp(sprintf('%s: skipping',sprintf('%s report %s',durationType{di},reportType{ci})));
            end
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
        info.call = 'nocv,bads';
        info.ci = ci;
        info.di = di;
        info.subj = -1; % all subjects
        info.dataType = sprintf('%s report %s',durationType{di},reportType{ci});
        infos{end+1} = info;
        disp(sprintf('ALLDATA %s: %i trials',info.dataType,size(data,1)));
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
    infos{ii}.fit = ac_fitTCCModel(infos{ii}.data,infos{ii}.call);
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
                        allfits{ci,di} = infos{ii}.fit;
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

%% Save the fits

save(fullfile('~/proj/afcom/tcc_data.mat'),'fits','allfits');

%% Load the fits

load(fullfile('~/proj/afcom/tcc_data.mat'));

%% Plot separated likelihood functions

for ci = 1
    for di = 2
        ac_plotSepLikes_tcc(allfits{ci,di});
    end
end

%% Load parameters from each subject
% pull out the dprime for spatial and feature (dt)
for ai = 1:5
    for ci = 1:2
        for di = 1:2
            fit = fits{ai}{ci,di};
            dt_all(ai,ci,di) = fit.params.dt_1;
            dt_s(ai,ci,di) = fit.params.dt_2;
            dt_f(ai,ci,di) = fit.params.dt_3;
            dt_baseline(ai,ci,di) = fit.params.dt_5;
        end
    end
    
end
% take the difference by condition
diff = dt_s-dt_f;
% reduce to the mean across difficulties and tasks for each subject
diff = squeeze(mean(mean(diff,3),2));
% bootstrap and plot
diff_ = mean(diff);
diff_ci = bootci(1000,@mean,diff);

% for plotting
dt_all = squeeze(mean(mean(dt_all,3),2));
dt_s = squeeze(mean(mean(dt_s,3),2));
dt_f = squeeze(mean(mean(dt_f,3),2));

dt_all_ci = bootci(1000,@mean,dt_all);
dt_s_ci = bootci(1000,@mean,dt_s);
dt_f_ci = bootci(1000,@mean,dt_f);

cmap = ac_cmap;
h = figure; hold on
errbar(1,mean(dt_all),dt_all_ci(2)-mean(dt_all),'-','Color',[0.75 0.75 0.75]);
plot(1,dt_all,'o','MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor','w');
plot(1,mean(dt_all),'o','MarkerFaceColor',[0.75 0.75 0.75],'MarkerEdgeColor','w','MarkerSize',10);
errbar(2,mean(dt_s),dt_s_ci(2)-mean(dt_s),'-','Color','k');
plot(2,dt_s,'o','MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor','w');
plot(2,mean(dt_s),'o','MarkerFaceColor','k','MarkerEdgeColor','w','MarkerSize',10);
errbar(3,mean(dt_f),dt_f_ci(2)-mean(dt_f),'-','Color','k');
plot(3,dt_f,'o','MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor','w');
plot(3,mean(dt_f),'o','MarkerFaceColor','k','MarkerEdgeColor','w','MarkerSize',10);
hline(mean(dt_baseline(:)),'--k');
axis([0 3 0 2]);
set(gca,'XTick',[1 2 3],'XTickLabel',{'Cue 4','Spatial','Feature'},'YTick',[0 1 2]);
ylabel('Sensitivity (d'')');
text(1,0.5,sprintf('Difference: %1.2f, 95%% CI [%1.2f, %1.2f]',diff_,diff_ci(1),diff_ci(2)));
drawPublishAxis('figSize=[20,10]','poster=1');


savepdf(h,fullfile('~/proj/afcom/figures','dprime.pdf'));