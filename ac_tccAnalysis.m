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
    
    %% Drop the first 100 trials
%     adata = adata(151:end,:);
%     disp(sprintf('Dropping 150 trials to account for training time'));

    
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
    mindurs = [0.25];
    maxdurs = [0.3];
    durationType = {'hard'};
    for ci = 1:length(cues)
        for di = 1:length(mindurs)
            data = sel(adata,3,cues(ci));
            data = fil(data,4,'>=',mindurs(di));
            data = fil(data,4,'<=',maxdurs(di));
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
    for di = 1:length(mindurs)
        data = sel(alldata,3,cues(ci));
        data = fil(data,4,'>=',mindurs(di));
        data = fil(data,4,'<=',maxdurs(di));
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
    infos{ii}.fit.data = infos{ii}.data;
    infos{ii}.fit.call = infos{ii}.call;
    infos{ii}.fit.dataType = infos{ii}.dataType;
end
disp('Fits complete');

%% Check that the fits 

%% Sort all the fits

allfits = {}; fits = {};
for si = [-1 1:length(subjects)]
    fit = {};
    for ci = 1:length(cues)
        for di = 1:length(mindurs)
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

save(fullfile('~/proj/afcom/tcc_data2.mat'),'fits','allfits');

%% Load the fits

load(fullfile('~/proj/afcom/tcc_data.mat'));

%% Plot the quality of the fits 
% using fits, a cell array of the fits to the two conditions, report
% color/direction

% Note that all of the non-target directions are uniformly distributed, so
% in practice you can show the fit of the model by just assuming that the
% target sensitivity is just reduced by 1-%target bias? We'll see if that
% works...

% average the respDistance for each condition separately (trialType

% TARGET TRIALTYPE CUE DURATION DEAD TARGETANGLE DISTRACTORANGLE
%    1      2       3     4       5       6           7
% ANGLE1 ANGLE2 ANGLE3 ANGLE4 RESPANGLE RESPDISTANCE DISTDISTANCE RUNS RT
%    8      9     10     11      12          13            14      15  16

% bin centers
xs = pi/64:pi/32:pi;

    reportType = {'color','direction'};
resp = zeros(length(fits),2,5,length(xs));
model = resp;
for subj = 1:length(fits)
    for cond = 1:2
        fit = fits{subj}{cond};

        for tt = 0:4
            ttdata = sel(fit.data,2,tt);
            % bin the values
            [c,~] = hist(ttdata(:,13),xs);
            % normalize counts (so we can average them later)
            c = c ./ sum(c);
            % save
            resp(subj,cond,tt+1,:) = c;
            
            % compute the model fit
            dt = fit.params.(sprintf('dt_%i',tt+1));
            probs = computeTCCPDF(xs,dt);
            model(subj,cond,tt+1,:) = probs;
        end
    end
end

%% Average everything and plot the fits

resp_mu = squeeze(nanmean(resp));
model_mu = squeeze(nanmean(model));

resp_ci = bootci(1000,@nanmean,resp);
model_ci = bootci(1000,@nanmean,model);

idxs = [1:12 14 16 20 24 32];

px = pscale(xs);

pos = [1 2 3 nan 4];
for cond = 1:2
    h = figure;
    for tt = [1 2 3 5]
        clear mmu mmu_ err err_ mu mu_
        subplot(4,1,pos(tt)); hold on
        px_ = px(idxs);
        mmu = squeeze(model_mu(cond,tt,:));
        mmu_ = mmu(idxs);
        plot(px_,mmu_,'-','Color',cmap(tt,:));
        err = squeeze(resp_ci(2,1,cond,tt,:))-squeeze(resp_mu(cond,tt,:));
        err_ = err(idxs);
        mu = squeeze(resp_mu(cond,tt,:));
        mu_ = mu(idxs);
        e = errbar(px_,mu_,err_,'-','Color',cmap(tt,:));
        plot(px_,mu_,'o','MarkerFaceColor',cmap(tt,:),'MarkerEdgeColor','w','MarkerSize',3);
        a = axis;
        axis([0 1 0 0.4]);
        set(gca,'XTick',[0 0.5 1],'YTick',[0 0.25]);
        if tt==5
            xlabel('Normalized psychphysical distance (a.u.)');
        ylabel('Response likelihood (pdf)');
        end
        drawPublishAxis('figSize=[4.4,8.9]');   
    end
%     savepdf(h,fullfile('~/proj/afcom/figures',sprintf('report%s_avg_model_fit.pdf',reportType{cond})));
end

%% Plot separated likelihood functions

for ci = 1:2
    for di = 1
        ac_plotSepLikes_tcc(allfits{ci,di});
    end
end

%% Load parameters from each subject
% pull out the dprime for spatial and feature (dt)

for ci_ = 1:2
    cmap_ = colorblindmap/255;

    cbmap(1,:) = cmap_(2,:);
    cbmap(2,:) = cmap_(3,:);

    % cmap_ = brighten(cmap_,.25);
    % 
    % cbmap(3,:) = cmap(2,:);
    % cbmap(4,:) = cmap(3,:);

    dt_all = zeros(8,2);
    bs_all = dt_all; bf_all = dt_all; bi_all = dt_all;
    dt_s = dt_all; bs_s = dt_all; bf_s = dt_all; bi_s = dt_all;
    for ai = 1:8
        for ci = 1:2
            fit = fits{ai}{ci};
            dt_all(ai,ci) = fit.params.dt_1;
            bs_all(ai,ci) = fit.params.bs_1;
            bf_all(ai,ci) = fit.params.bf_1;
            bi_all(ai,ci) = fit.params.bi_1;
            dt_s(ai,ci) = fit.params.dt_2;
            bs_s(ai,ci) = fit.params.bs_2;
            bf_s(ai,ci) = fit.params.bf_2;
            bi_s(ai,ci) = fit.params.bi_2;
            dt_f(ai,ci) = fit.params.dt_3;
            bs_f(ai,ci) = fit.params.bs_3;
            bf_f(ai,ci) = fit.params.bf_3;
            bi_f(ai,ci) = fit.params.bi_3;
            dt_baseline(ai,ci) = fit.params.dt_5;
        end
    end

% normalize dt_s and dt_f to the baseline
% dt_all = dt_all ./ dt_baseline;
% dt_s = dt_s ./ dt_baseline;
% dt_f = dt_f ./ dt_baseline;

% Plot this twice -- once for the report color task, once for the report direction task

    dt_all = dt_all(:,ci_);
    dt_s = dt_s(:,ci_);
    dt_f = dt_f(:,ci_);
    dt_baseline = dt_baseline(:,ci_);
    
    dt_all = dt_all ./ dt_baseline;
    dt_s = dt_s ./ dt_baseline;
    dt_f = dt_f ./ dt_baseline;

    % take the difference by condition
    diff = dt_s-dt_f;
    % reduce to the mean across difficulties and tasks for each subject
    diff = squeeze(mean(mean(diff,3),2));
    % bootstrap and plot
    diff_ = mean(diff);
    diff_ci = bootci(1000,@mean,diff);

    % for plotting
    dt_all_ = squeeze(mean(mean(dt_all,3),2));
    dt_s_ = squeeze(mean(mean(dt_s,3),2));
    dt_f_ = squeeze(mean(mean(dt_f,3),2));

    dt_all_ci = bootci(1000,@mean,dt_all);
    dt_s_ci = bootci(1000,@mean,dt_s);
    dt_f_ci = bootci(1000,@mean,dt_f);


    cmap = ac_cmap;
    h = figure; hold on
    
    % plot individuals as lines
    for ai = 1:size(dt_all,1)
        plot([1 2 3],[dt_all(ai) dt_s(ai) dt_f(ai)],'-','Color',[0.75 0.75 0.75]);
    end
    
    errbar(1,mean(dt_all_),dt_all_ci(2)-mean(dt_all_),'-','Color',[0.75 0.75 0.75]);
    plot(1,dt_all_,'o','MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor','w','MarkerSize',5);
    plot(1,mean(dt_all_),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor','w','MarkerSize',7);
    errbar(2,mean(dt_s_),dt_s_ci(2)-mean(dt_s_),'-','Color','k');

    % spatial
    plot(2,dt_s_,'o','MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor','w','MarkerSize',5);
    plot(2,mean(dt_s_),'o','MarkerFaceColor','k','MarkerEdgeColor','w','MarkerSize',7);

    % feature
    errbar(3,mean(dt_f_),dt_f_ci(2)-mean(dt_f_),'-','Color','k');
    plot(3,dt_f_,'o','MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor','w','MarkerSize',5);
    plot(3,mean(dt_f_),'o','MarkerFaceColor','k','MarkerEdgeColor','w','MarkerSize',7);
    axis([0 3 0 1]);
    set(gca,'XTick',[1 2 3],'XTickLabel',{'Cue 4','Cue side','Cue feature'},'YTick',[0 1]);
    ylabel('Sensitivity (d'')');
    text(0,1,'Ceiling');
    hline(1,'--k');
%     hline(mean(dt_baseline(:)),'--k');
%     text(0,mean(dt_baseline(:)),'Ceiling');
    text(1,0.5,sprintf('Difference: %1.2f, 95%% CI [%1.2f, %1.2f]',diff_,diff_ci(1),diff_ci(2)));
    drawPublishAxis('figSize=[4.5,4.5]');


% 	savepdf(h,fullfile('~/proj/afcom/figures',sprintf('dprime_%s.pdf',reportType{ci_})));
end

%% Save data again
dt_all = zeros(8,2);
bs_all = dt_all; bf_all = dt_all; bi_all = dt_all;
dt_s = dt_all; bs_s = dt_all; bf_s = dt_all; bi_s = dt_all;
for ai = 1:8
    for ci = 1:2
        fit = fits{ai}{ci};
        dt_all(ai,ci) = fit.params.dt_1;
        bs_all(ai,ci) = fit.params.bs_1;
        bf_all(ai,ci) = fit.params.bf_1;
        bi_all(ai,ci) = fit.params.bi_1;
        dt_s(ai,ci) = fit.params.dt_2;
        bs_s(ai,ci) = fit.params.bs_2;
        bf_s(ai,ci) = fit.params.bf_2;
        bi_s(ai,ci) = fit.params.bi_2;
        dt_f(ai,ci) = fit.params.dt_3;
        bs_f(ai,ci) = fit.params.bs_3;
        bf_f(ai,ci) = fit.params.bf_3;
        bi_f(ai,ci) = fit.params.bi_3;
        dt_baseline(ai,ci) = fit.params.dt_5;
    end
end
%% Also generate a plot which shows the fit of the TCC model (w/ these parameters)
% to the individual subject data? 

%% Now generate plots showing how bias changes between conditions

to = -0.025;
sgap = 0.005;
bgap = 0.01;

for cond = 1:2
    h = figure;

    % generate three subplots
    subplot(3,1,1); hold on
    
    % ALL CONDITION
    
    % average the values
    bs = nanmean(bs_all(:,cond));
    bf = nanmean(bf_all(:,cond));
    bi = nanmean(bi_all(:,cond));
    
    subplot(311); hold on
    % this helper function assumes there is a value "vals" that exists
    % which has the bs/bf/bi beta weights in it. It does the rest of the
    % job of plotting them into the axes with the appropriate scale
    help_ac_plot_bias;
    
    bs = nanmean(bs_s(:,cond));
    bf = nanmean(bf_s(:,cond));
    bi = nanmean(bi_s(:,cond));
    subplot(312); hold on
    help_ac_plot_bias;
    
    
    bs = nanmean(bs_f(:,cond));
    bf = nanmean(bf_f(:,cond));
    bi = nanmean(bi_f(:,cond));
    subplot(313); hold on
    help_ac_plot_bias;
    
    reportType = {'color','direction'};
%     savepdf(h,fullfile('~/proj/afcom/figures',sprintf('ac_bias_%s.pdf',reportType{cond})));
    
    % TODO: Output based on condition, and add labels to the different
    % conditions (or do this in the paper? 
end


%% Now generate a bias plot that is the average of the two conditions
to = -0.025;
sgap = 0.005;
bgap = 0.01;

h = figure;

% generate three subplots
subplot(3,1,1); hold on

% ALL CONDITION

% average the values
bs = nanmean(bs_all(:));
bf = nanmean(bf_all(:));
bi = nanmean(bi_all(:));

subplot(311); hold on
% this helper function assumes there is a value "vals" that exists
% which has the bs/bf/bi beta weights in it. It does the rest of the
% job of plotting them into the axes with the appropriate scale
help_ac_plot_bias;

bs = nanmean(bs_s(:));
bf = nanmean(bf_s(:));
bi = nanmean(bi_s(:));
subplot(312); hold on
help_ac_plot_bias;


bs = nanmean(bs_f(:));
bf = nanmean(bf_f(:));
bi = nanmean(bi_f(:));
subplot(313); hold on
help_ac_plot_bias;

reportType = {'color','direction'};
savepdf(h,fullfile('~/proj/afcom/figures',sprintf('ac_bias_average.pdf')));
    


%% Statistics

%% Save data again
dt_all = zeros(8,2);
bs_all = dt_all; bf_all = dt_all; bi_all = dt_all;
dt_s = dt_all; bs_s = dt_all; bf_s = dt_all; bi_s = dt_all;
for ai = 1:8
    for ci = 1:2
        fit = fits{ai}{ci};
        dt_all(ai,ci) = fit.params.dt_1;
        bs_all(ai,ci) = fit.params.bs_1;
        bf_all(ai,ci) = fit.params.bf_1;
        bi_all(ai,ci) = fit.params.bi_1;
        dt_s(ai,ci) = fit.params.dt_2;
        bs_s(ai,ci) = fit.params.bs_2;
        bf_s(ai,ci) = fit.params.bf_2;
        bi_s(ai,ci) = fit.params.bi_2;
        dt_f(ai,ci) = fit.params.dt_3;
        bs_f(ai,ci) = fit.params.bs_3;
        bf_f(ai,ci) = fit.params.bf_3;
        bi_f(ai,ci) = fit.params.bi_3;
        dt_baseline(ai,ci) = fit.params.dt_5;
    end
end

%%

% we've pulled out all the parameters here, now we need to average them and
% compare them across conditions 

% for example, is there a difference between dt for cue side and cue
% feature? 
% start by computing on aggregate data
diff = dt_f-dt_s;
ci = bootci(10000,@nanmean,diff(:));
disp(sprintf('Cue Feature was %1.2f $d''$ higher, 95\\%% CI [%1.2f, %1.2f]',mean(diff(:)),ci(1),ci(2)));

