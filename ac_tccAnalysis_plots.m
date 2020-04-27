
%% Load the fits

load(fullfile('~/proj/afcom/tcc_data.mat'));

%% Setup colors
% we'll use colorblind colors, which has 8 colors
% we need four for the model comparison (2-4 and 8)

% we need 3 for the dprime comparison (5-7)

cmap = colorblindmap/255;
% set the colors for the task
colors = struct;
colors.all = cmap(5,:);
colors.baseline = cmap(6,:);
colors.shared = cmap(7,:);
colors.spatial = cmap(2,:);
colors.feature = cmap(3,:);

%% Remove data that wasn't actually used in fitting
trialTypes = {'all','spatial','feature','target','baseline'};
trialTypeVals = [0 1 2 3 4];

for si = 1:8
    for ci = 1:2
        for mi = 1:6
            fit = fits{si}{ci,mi};
            adata = fit.data;
            
            cTrialTypes = trialTypes;
            cTrialTypeVals = trialTypeVals;
            
            mode = fit.call;
            
            keep = cellfun(@(x) ~isempty(strfind(mode,x)),cTrialTypes);
            cTrialTypes = cTrialTypes(keep);
            cTrialTypeVals = cTrialTypeVals(keep);

            % remove trials we aren't fitting
            if length(cTrialTypes)<5
                disp(sprintf('Removing %i conditions that are not being fit',5-length(cTrialTypes)));
                cTrialTypes = adata(:,2);
                idxs = zeros(size(cTrialTypes));
                for tt = 1:length(cTrialTypeVals)
                    idxs = idxs + (cTrialTypes==cTrialTypeVals(tt));
                end
                disp(sprintf('Dropped %i trials',size(adata,1)-sum(idxs>0)));
                adata = adata(logical(idxs),:);
            end
            
            fit.data = adata;
            clear adata
            fits{si}{ci,mi} = fit;
        end
    end
end

%% Get the permutations and look at their distributions
perms = nan(8,2,6,100);
like = nan(8,2,6);

for si = 1:length(fits)
    for ci = 1:2
        for mi = 1:6
            cfit = fits{si}{ci,mi};
            if isfield(cfit,'perm') && all(cfit.perm~=0)
                perms(si,ci,mi,:) = cfit.perm;
            end
            like(si,ci,mi) = cfit.cv.likelihood;
        end
    end
end

%% take all the permutations and remove the mean for each set
perms_ = perms - repmat(nanmean(perms,4),1,1,1,100);
% now find the 95% quantiles
ci_perms = quantile(perms_(:),[.025 .975]);

%% Check that the model likelihoods are > permutation likelihoods
perm_mean = nanmean(perms,4);
diff = perm_mean-like;

figure;
hist(diff(:));

% indeed, all models except one are better, and that one is equal. I.e.
% probably the data is garbage from that participant and they couldn't do
% that task?

%%
h = figure;

% just plot them all side by side
num = [];
for si = 1:size(perms,1)
    for ci = 1:size(perms,2)
        for mi = 1:size(perms,3)
            vals = perms(si,ci,mi,:);
            if any(~isnan(vals))
                ci95 = bootci(1000,@nanmean,vals);
                num(1,end+1) = ci95(1);
                num(2,end) = ci95(2);
            end
        end
    end
end

% the point of this is to show that the spread of permutations is tiny (<1
% in the cross-validated information criterion). This means that any
% differences >1 are probably significant, so the model comparisons we've
% done are almost all useful to look at. It also means that the actual
% models which are all WAY better fits than the permutations were valuable.
hist(num(2,:)-num(1,:));


%% Cross-validated likelihood
% Let's pull all the CV likelihoods out to compare them, specifically the
% cross-validated fits for the 4 models
cv = [];
for si = 1:length(fits)
    for ci = 1:2
        for mi = 3:6
            cv(si,ci,mi-2) = fits{si}{ci,mi}.cv.likelihood;
        end
    end
end

% now we can check whether the numbers are bigger or smaller. If a model
% has a *smaller* value, itis a better model

% CV model #:
%  1         2    3      4
% all sh   bias sens bias+sens

% we're going to check if 

 % critical comparison: if greater than zero, then there is value to using
 % non-shared bias parameters
temp1 = cv(:,1,1)-cv(:,1,2);
bootci(10000,@mean,temp1(:))
% also do a ranksum test

% if > 0, then there is value to using non-shared sensitivity parameters
temp2 = cv(:,1,1)-cv(:,1,3);
bootci(10000,@mean,temp2(:))

% if > 0 then there is value to using non-shared parameters
temp3 = cv(:,1,1)-cv(:,1,4);
bootci(10000,@mean,temp3(:))

%% Signed rank tests for the cross-validated likelihoods

% check whether each set of CVs is better than the baseline model, this
% means that the cv(:,1,1) is LARGER than the better model, then it's
% better (or whatever)
cv_demeaned = repmat(cv(:,:,1),1,1,4) - cv;

for cond = 1:2
    for compare = 2:4
        p = signrank(squeeze(cv_demeaned(:,cond,compare)),0,'tail','right');
        disp(p);
    end
end


%% Mean cross-validated likelihoods
% remove the means
cv_demeaned = repmat(cv(:,:,1),1,1,4) - cv;

% get the mean across subjects of each condition
cv_ = squeeze(mean(cv_demeaned,1));
cv_ci95 = squeeze(bootci(1000,@nanmean,cv_demeaned));
cv_ci95 = squeeze(cv_ci95(2,:,:)) - cv_;

h = figure;

titles = {'Selection direction, report color','Select color, report direction'};

for cond = 1:2
    subplot(2,1,cond); hold on
    plot([2 4],[0 0],'--k'); %,'Color',cmap(1,:));
    for mi = 2:4
        % add individuals first
        plot(mi,cv_demeaned(:,cond,mi),'o','MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor','w','MarkerSize',5);
        
        errbar(mi,cv_(cond,mi),cv_ci95(cond,mi),'-','Color','k');
        plot(mi,cv_(cond,mi),'o','MarkerFaceColor','k','MarkerEdgeColor','w','MarkerSize',8);
        % color version
%         errbar(mi,cv_(cond,mi),cv_ci95(cond,mi),'-','Color',cmap(mi,:));
%         plot(mi,cv_(cond,mi),'o','MarkerFaceColor',cmap(mi,:),'MarkerEdgeColor','w');
    end

    % permutation test significance
%     temp = hline(ci_perms(1),'--r');
%     set(temp,'Color',cmap(8,:));
%     temp = hline(ci_perms(2),'--r');
%     set(temp,'Color',cmap(8,:));

    % just use +/- 2
    hline(-2,'--r');
    hline(2,'--r');

    title(titles{cond});
    axis([2 4 -10 15]);
    
    %     xlabel('Observer');
    ylabel('Relative fit (\Delta cross-validated likelihood)');

    %     legend(p,{'Shared parameters','Different bias','Different sensitivity','No shared parameters'});

    set(gca,'XTick',[2 3 4],'XTickLabel',{'Separate bias','Separate sens','Both separate'},'YTick',[-5 0 5],'YTickLabel',{'Worse -5','0','Better +5'});

    drawPublishAxis('figSize=[8.9,10]');

end

savepdf(h,fullfile('~/proj/afcom/figures','modelComparison_avg.pdf'));

%% Plot the overall cross-validated likelihoods for the different models in 
% a way that you can compare for each subject individually

h = figure;

subplot(211); hold on

cmap = colorblindmap/255;

title('Select direction, report color');

offsets = [1 2 3 4]; offsets = offsets - mean(offsets); offsets = offsets/8;

p(1) = plot([1 8],[0 0],'--','Color',cmap(1,:));

% add two more hlines which are the 95% confidence intervals from the
% permutation tests. In fact since these are ~2, we'll just use 2 and be
% consistent with the literature.
temp = hline(-2,'--');
set(temp,'Color',cmap(8,:));
temp = hline(2,'--');
set(temp,'Color',cmap(8,:));

for si = 1:8
    for mi = 2:4
        like = cv(si,1,1) - cv(si,1,mi);
        
        p(mi) = plot(si+offsets(mi-1),like,'o','MarkerFaceColor',cmap(mi,:),'MarkerEdgeColor','w');
    end
end

xlabel('Observer');
ylabel('Relative fit (\Delta cross-validated likelihood)');

legend(p,{'Shared parameters','Different bias','Different sensitivity','No shared parameters'});

set(gca,'XTick',1:8,'YTick',[-5 0 5 10],'YTickLabel',{'Worse -5','0','+5','Better +10'});

drawPublishAxis;

subplot(212); hold on
title('Select color, report direction');

p(1) = plot([1 8],[0 0],'--','Color',cmap(1,:));

% add two more hlines which are the 95% confidence intervals from the
% permutation tests
temp = hline(ci_perms(1),'--');
set(temp,'Color',cmap(8,:));
temp = hline(ci_perms(2),'--');
set(temp,'Color',cmap(8,:));

for si = 1:8
    for mi = 2:4
        like = cv(si,2,1) - cv(si,2,mi);
        
        p(mi) = plot(si+offsets(mi-1),like,'o','MarkerFaceColor',cmap(mi,:),'MarkerEdgeColor','w');
    end
end

xlabel('Observer');
ylabel('Relative fit (\Delta cross-validated likelihood)');

legend(p,{'Shared parameters','Different bias','Different sensitivity','No shared parameters'});

set(gca,'XTick',1:8,'YTick',[-5 0 5 10],'YTickLabel',{'Worse -5','0','+5','Better +10'});

drawPublishAxis('figSize=[15,8.9]');

savepdf(h,fullfile('~/proj/afcom/figures/','modelComparison.pdf'));

%% Plot the quality of the fits 
% using fits, a cell array of the fits to the two conditions, report
% color/direction

% Note that all of the non-target directions are uniformly distributed, so
% in practice you can show the fit of the model by just assuming that the
% target sensitivity is just reduced by 1-%target bias? We'll see if that
% works... <- note from 3/20/20, I have no idea what that comment is
% about?

% average the respDistance for each condition separately (trialType

% TARGET TRIALTYPE CUE DURATION DEAD TARGETANGLE DISTRACTORANGLE
%    1      2       3     4       5       6           7
% ANGLE1 ANGLE2 ANGLE3 ANGLE4 RESPANGLE RESPDISTANCE DISTDISTANCE RUNS RT
%    8      9     10     11      12          13            14      15  16

% Note that because the wrong data was kept in the fit, we have to remove
% data types that are irrelevant

% bin centers
xs = pi/64:pi/32:pi;
% actual bins
xbins = 0:pi/32:pi;

reportType = {'color','direction'};
resp = nan(length(fits),2,6,5,length(xs));
model = resp;
for subj = 1:length(fits)
    for cond = 1:2
        for mi = 1:6
            % pull the current fit object
            cfit = fits{subj}{cond,mi};
            
            % for each trial type, bin the response distances
            tts = unique(cfit.data(:,2));
            for ti = 1:length(tts)
                tt = tts(ti);
                
                % bin the values
                tdata = cfit.data(cfit.data(:,2)==tt,:);
                
                [c] = histc(tdata(:,13),xbins);
                c(end-1) = c(end-1)+c(end); % in case any values exactly match pi
                c = c(1:end-1); % remove the last value which is now empty
                
                % normalize data for averaging later
                c = c ./ sum(c);
                % save
                resp(subj,cond,mi,tt+1,:) = c;
                
                % get the model fits, note that the choice of dt depends on
                % the type of model that was fit here
                if isfield(cfit.params,sprintf('dt_%i',tt))
                    dt = cfit.params.(sprintf('dt_%i',tt));
                else
                    dt = cfit.params.dt_sh;
                end
                probs = computeTCCPDF(xs,dt);
                if any(isnan(probs))
                    disp('wtf')
                end
                model(subj,cond,mi,tt+1,:) = probs;
            end
        end
    end
end

%% Let's compute R^2 on these binned responses while we're here

% there are two conditions and there are six models, compute r^2 for each
% one

% For r^2 we will use the 1-res/total approach, appropriate for a nonlinear
% fit
SSres = sum((resp-model).^2,5); % mean squared error
SStotal = sum(resp.^2,5);

r2 = 1-SSres./SStotal;
% r2 = 1 - (SSres./SStotal);

% correlation option for r^2
% r2 = zeros(size(resp,1),size(resp,2),size(resp,3),size(resp,4));
% for i = 1:8
%     for j = 1:2
%         for k = 1:6
%             for l = 1:6
%                 r2(i,j,j,l) = corr(squeeze(model(i,j,k,l,:)),squeeze(resp(i,j,k,l,:))).^2;
%             end
%         end
%     end
% end

r2 = r2(:,:,:,[1 2 3 5]); % remove trial types that aren't uncued, spatial, feature, or baseline

% get the r^2 for the all condition and the baseline condition 
r2_uncued = squeeze(r2(:,:,1,1));
for cond = 1:2
    mu = mean(r2_uncued(:,cond));
    ci = bootci(1000,@nanmean,r2_uncued(:,cond));
    disp(sprintf('%s: %1.2f [%1.2f, %1.2f]',reportType{cond},mu,ci(1),ci(2)));
end
r2_nodist = squeeze(r2(:,:,2,4));
for cond = 1:2
    mu = mean(r2_nodist(:,cond));
    ci = bootci(1000,@nanmean,r2_nodist(:,cond));
    disp(sprintf('%s: %1.2f [%1.2f, %1.2f]',reportType{cond},mu,ci(1),ci(2)));
end
disp('spatial');
r2_spatial = squeeze(r2(:,:,3,2));
for cond = 1:2
    mu = mean(r2_spatial(:,cond));
    ci = bootci(1000,@nanmean,r2_spatial(:,cond));
    disp(sprintf('%s: %1.2f [%1.2f, %1.2f]',reportType{cond},mu,ci(1),ci(2)));
end
disp('feature');
r2_feature = squeeze(r2(:,:,3,3));
for cond = 1:2
    mu = mean(r2_feature(:,cond));
    ci = bootci(1000,@nanmean,r2_feature(:,cond));
    disp(sprintf('%s: %1.2f [%1.2f, %1.2f]',reportType{cond},mu,ci(1),ci(2)));
end

%% Now compare the R^2 improvement that results from forcing shared bias vs. shared sensitivity
% i.e. compute the full separate r^2 minus shared bias or minus shared
% sensitivity. These changes are reported in the paper for the last figure
% as a sort of "effect size" calculation. 

%% compare the weird ones (subj 5)
s5resp = squeeze(resp(5,1,1,1,:));
s5model = squeeze(model(5,1,1,1,:));

h =figure; hold on
plot(xs,s5resp,'o');
plot(xs,s5model,'-');

%% Average everything and plot the fits
cmap = colorblindmap/255;

cmap(4,:) = [0 0 0];

% use model #3 the shared fit

resp_use(:,:,1,:) = resp(:,:,1,1,:);
resp_use(:,:,2,:) = resp(:,:,3,2,:);
resp_use(:,:,3,:) = resp(:,:,4,3,:);
resp_use(:,:,4,:) = resp(:,:,2,5,:);

model_use(:,:,1,:) = model(:,:,1,1,:);
model_use(:,:,2,:) = model(:,:,3,2,:);
model_use(:,:,3,:) = model(:,:,4,3,:);
model_use(:,:,4,:) = model(:,:,2,5,:);

resp_mu = squeeze(nanmean(resp_use));
model_mu = squeeze(nanmean(model_use));

resp_ci = bootci(1000,@nanmean,resp_use);
model_ci = bootci(1000,@nanmean,model_use);

idxs = [1:12 14 16 20 24 32];

px = pscale(xs);

pos = [1 2 3 4];
tts = {0 1 2 4};
titles = {'Cue 4','Cue spatial','Cue feature','Baseline'};
for cond = 1:2
    h = figure;
    for ti = [1 2 3 4]
        tt = tts{ti};
        clear mmu mmu_ err err_ mu mu_
        subplot(max(pos),1,pos(ti)); hold on
        title(titles{pos(ti)});
        px_ = px(idxs);
        mmu = squeeze(model_mu(cond,ti,:));
        mmu_ = mmu(idxs);
        plot(px_,mmu_,'-','Color',cmap(ti,:));
        err = squeeze(resp_ci(2,1,cond,ti,:))-squeeze(resp_mu(cond,ti,:));
        err_ = err(idxs);
        mu = squeeze(resp_mu(cond,ti,:));
        mu_ = mu(idxs);
        e = errbar(px_,mu_,err_,'-','Color',cmap(ti,:));
        plot(px_,mu_,'o','MarkerFaceColor',cmap(ti,:),'MarkerEdgeColor','w','MarkerSize',3);
        a = axis;
        axis([0 1 0 0.4]);
        set(gca,'XTick',[0 0.5 1],'YTick',[0 0.25]);
        if tt==4
            xlabel('Normalized psychphysical distance (a.u.)');
        ylabel('Response likelihood (pdf)');
        end
        drawPublishAxis('figSize=[4.4,8.9]');   
    end
    savepdf(h,fullfile('~/proj/afcom/figures',sprintf('report%s_avg_model_fit.pdf',reportType{cond})));
end

%% Get the dprime parameters and plot these against the baseline/all conditions

dps = zeros(8,2,6);
param_name = {'dt_sh','dt_sh','dt_sh','dt_sh'};

% the dp parameter can be shared or specific, separate those
dps_ab = zeros(8,2,4);
dps_sf = zeros(8,2,2,2);

for subj = 1:length(fits)
    for cond = 1:2
        for mi = 1:4
            cfit = fits{subj}{cond,mi};
            
            dps_ab(subj,cond,mi) = cfit.params.dt_sh;
        end
        
        % models where the dp is NOT shared 
        for mi = 5:6
            cfit = fits{subj}{cond,mi};
            
            dps_sf(subj,cond,mi-4,1) = cfit.params.dt_1;
            dps_sf(subj,cond,mi-4,2) = cfit.params.dt_2;
        end
    end
end

%% get the mean and CI for various parts of the paper
mu = mean(dps_ab(:,1,1));
ci = bootci(1000,@mean,dps_ab(:,1,1));
disp(sprintf('Uncued mean report-color %1.2f 95%% CI [%1.2f, %1.2f]',mu,ci(1),ci(2)));
mu = mean(dps_ab(:,2,1));
ci = bootci(1000,@mean,dps_ab(:,2,1));
disp(sprintf('Uncued mean report-direction %1.2f 95%% CI [%1.2f, %1.2f]',mu,ci(1),ci(2)));

mu = mean(dps_ab(:,1,2));
ci = bootci(1000,@mean,dps_ab(:,1,2));
disp(sprintf('Baseline mean report-color %1.2f 95%% CI [%1.2f, %1.2f]',mu,ci(1),ci(2)));
mu = mean(dps_ab(:,2,2));
ci = bootci(1000,@mean,dps_ab(:,2,2));
disp(sprintf('Baseline mean report-direction %1.2f 95%% CI [%1.2f, %1.2f]',mu,ci(1),ci(2)));

% test baseline against uncued
diff = squeeze(dps_ab(:,1,2)-dps_ab(:,1,1));
mu = mean(diff);
ci = bootci(10000,@nanmean,diff);
disp(sprintf('No-distractor - uncued report-color mu = %1.2f, 95%% CI [%1.2f, %1.2f]',mu,ci(1),ci(2)));
diff = squeeze(dps_ab(:,2,2)-dps_ab(:,2,1));
mu = mean(diff);
ci = bootci(10000,@nanmean,diff);
disp(sprintf('No-distractor - uncued report-direction mu = %1.2f, 95%% CI [%1.2f, %1.2f]',mu,ci(1),ci(2)));

mu = mean(dps_ab(:,1,3));
ci = bootci(1000,@mean,dps_ab(:,1,3));
disp(sprintf('Shared mean report-color %1.2f, 95%% CI [%1.2f, %1.2f]',mu,ci(1),ci(2)));
mu = mean(dps_ab(:,2,3));
ci = bootci(1000,@mean,dps_ab(:,2,3));
disp(sprintf('Shared mean report-direction %1.2f, 95%% CI [%1.2f, %1.2f]',mu,ci(1),ci(2)));

% difference between shared mean and uncued mean
shmu = dps_ab(:,:,3);
unmu = dps_ab(:,:,1);

% tail right checks if X > Y
diff = squeeze(shmu(:,1)-unmu(:,1));%,0,'tail','right');
mu = mean(diff);
ci = bootci(10000,@nanmean,diff);
disp(sprintf('Shared - uncued report-color mu = %1.2f, 95%% CI [%1.2f, %1.2f]',mu,ci(1),ci(2)));
diff = squeeze(shmu(:,2)-unmu(:,2));%,0,'tail','right');
mu = mean(diff);
ci = bootci(10000,@nanmean,diff);
disp(sprintf('Shared - uncued report-direction mu = %1.2f, 95%% CI [%1.2f, %1.2f]',mu,ci(1),ci(2)));

% check if separate s/f parameters are different from the shared parameter
spatmu = dps_sf(:,:,2,1);
featmu = dps_sf(:,:,2,2);
[p,h,stats] = signrank(spatmu(:,1)-shmu(:,1));
disp(sprintf('Spatial report-color unshared - shared signed rank P=%1.2f',p));
[p,h,stats] = signrank(spatmu(:,2)-shmu(:,2));
disp(sprintf('Spatial report-direction unshared - shared signed rank P=%1.2f',p));
[p,h,stats] = signrank(featmu(:,1)-shmu(:,1));
disp(sprintf('Feature report-color unshared - shared signed rank P=%1.2f',p));
[p,h,stats] = signrank(featmu(:,2)-shmu(:,2));
disp(sprintf('Feature report-direction unshared - shared signed rank P=%1.2f',p));

% normalize dprime)
normalized = max(zeros(size(dps_ab(:,:,1))),(dps_ab(:,:,3)-dps_ab(:,:,1)))./(dps_ab(:,:,2)-dps_ab(:,:,1));

mu = mean(normalized);
ci = bootci(10000,@nanmean,normalized);
disp(sprintf('Report-color normalized dp %1.1f%%, 95%% CI [%1.1f, %1.1f]',mu(1)*100,ci(1,1)*100,ci(2,1)*100));
disp(sprintf('Report-direction normalized dp %1.1f%%, 95%% CI [%1.1f, %1.1f]',mu(2)*100,ci(1,2)*100,ci(2,2)*100));
%% Plot the dp of the shared fit (which is basically the best) 
% against the baseline and all fits, to get a sense of how much improvement
% happens


h = figure;

reportType = {'color','direction'};

% take the mean across subjs
dps_ab_ = squeeze(mean(dps_ab,1));
dps_ab_ci = squeeze(bootci(1000,@nanmean,dps_ab));

dps_sf_ = squeeze(mean(dps_sf,1));
dps_sf_ci = squeeze(bootci(1000,@nanmean,dps_sf));

% plot the color fits
for cond = 1:2
    subplot(2,1,cond); hold on
    title(sprintf('Report %s',reportType{cond}));

    color_opts = {colors.all colors.baseline colors.shared};
    color_opts2 = {colors.spatial colors.feature};
    pos = [1 5 2];
    pos2 = [3 4];

    for subj = 1:8
%         plot(pos,squeeze(dps_ab(subj,cond,1:3)),'-','Color',[0.5 0.5 0.5]);
        for i = 1:3
            plot(pos(i),dps_ab(subj,cond,i),'o','MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor','w');
        end
        for j = 1:2
            % the 2 in the third column is to grab the "full" model with
            % all parameters fit
            plot(pos2(j),dps_sf(subj,cond,2,j),'o','MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor','w');
        end
    end

    clear p
    for i = 1:3
        p(pos(i)) = plot(pos(i),dps_ab_(cond,i),'o','MarkerFaceColor',color_opts{i},'MarkerEdgeColor','w'); % all
    end
    for j = 1:2
        p(pos2(j)) = plot(pos2(j),dps_sf_(cond,2,j),'o','MarkerFaceColor',color_opts2{j},'MarkerEdgeColor','w'); % all
    end

%     legend(p,{'All','Shared','Spatial','Feature','Baseline'},'Location','NorthWest');

    axis([-1 7 0 3]);

    set(gca,'XTick',[1 2 3 4 5],'XTickLabel',{'All','Shared','Spatial','Feature','Baseline'});

    ylabel('Sensitivity (d'')');

    drawPublishAxis('figSize=[4.5 8.9]');
end

savepdf(h,fullfile('~/proj/afcom/figures/dp_comparison.pdf'));


%% Get the bias parameters and plot these against the baseline/all conditions

dps = zeros(8,2,6);

% the dp parameter can be shared or specific, separate those
bs = nan(8,2,3);
bf = nan(8,2,3);
bi = nan(8,2,3);
% specific parameters (for the spatial/feature fits)
bs_sf = nan(8,2,2,2);
bf_sf = nan(8,2,2,2);
bi_sf = nan(8,2,2,2);

for subj = 1:length(fits)
    for cond = 1:2
        % models with shared bias parameters
        idx = [1 2 3 -1 4];
        for mi = [1 2 3 5]
            cfit = fits{subj}{cond,mi};
            
            bs(subj,cond,idx(mi)) = cfit.params.bs_sh;
            bf(subj,cond,idx(mi)) = cfit.params.bf_sh;
            bi(subj,cond,idx(mi)) = cfit.params.bi_sh;
        end
        
        % models where the bias params are NOT shared 
        idx = [-1 -1 -1 1 -1 2];
        for mi = [4 6]
            cfit = fits{subj}{cond,mi};
            
            bs_sf(subj,cond,idx(mi),1) = cfit.params.bs_1;
            bs_sf(subj,cond,idx(mi),2) = cfit.params.bs_2;
            
            bf_sf(subj,cond,idx(mi),1) = cfit.params.bf_1;
            bf_sf(subj,cond,idx(mi),2) = cfit.params.bf_2;
            
            bi_sf(subj,cond,idx(mi),1) = cfit.params.bi_1;
            bi_sf(subj,cond,idx(mi),2) = cfit.params.bi_2;
        end
    end
end

%% Pull the values directly. The plots are kind of shit to work with
% we'll just build the plots in illustrator directly instead

% take the mean across subjects for the fits

bs_ = squeeze(mean(bs));
bf_ = squeeze(mean(bf));
bi_ = squeeze(mean(bi));

models = {'all','baseline','shared','shared_bias'};
reportType = {'color','direction'};

for cond = 1:2
    
    for mi = 1:4
        
        cbs = bs_(cond,mi);
        cbf = bf_(cond,mi);
        cbi = bi_(cond,mi);
        
        target = cbs * cbf;
        side = cbs * (1-cbf);
        feature = (1-cbs) * (1-cbi);
        distractor = (1-cbs) * cbi;
        
        disp(sprintf('Report %s, model %s',reportType{cond},models{mi}));
        disp(sprintf('Target %1.2f%%',target*100));
        disp(sprintf('Side %1.2f%%',side*100));
        disp(sprintf('Feature %1.2f%%',feature*100));
        disp(sprintf('Distractor %1.2f%%',distractor*100));
    end
end

%% Plot the bias on comparison figures specific for each parameter across models, like the sensitivity
% we first need to convert the bs/bf/bi into target/side/feature/distractor
% we also only need to get five values:
% uncued (cue 4)
% shared
% spatial (from full model)
% feature (from full model)
% no-distractor (all)
% columns are: subj | cond | model | bias parameter
bias_data = nan(8,2,5,4);

% make some temp holders
tt = bs .* bf;
ts = bs .* (1-bf);
tf = (1-bs) .* (1-bi);
ti = (1-bs) .* bi;

% copy values into full data
% uncued
bias_data(:,:,1,1) = tt(:,:,1);
bias_data(:,:,1,2) = ts(:,:,1);
bias_data(:,:,1,3) = tf(:,:,1);
bias_data(:,:,1,4) = ti(:,:,1);

% no-distractor / all
bias_data(:,:,5,1) = tt(:,:,2);
bias_data(:,:,5,2) = ts(:,:,2);
bias_data(:,:,5,3) = tf(:,:,2);
bias_data(:,:,5,4) = ti(:,:,2);

% shared
bias_data(:,:,2,1) = tt(:,:,3);
bias_data(:,:,2,2) = ts(:,:,3);
bias_data(:,:,2,3) = tf(:,:,3);
bias_data(:,:,2,4) = ti(:,:,3);

% now use the spatial/feature values
tt_sf = bs_sf .* bf_sf;
ts_sf = bs_sf .* (1-bf_sf);
tf_sf = (1-bs_sf) .* (1-bi_sf);
ti_sf = (1-bs_sf) .* bi_sf;

% spatial
bias_data(:,:,3,1) = tt_sf(:,:,2,1);
bias_data(:,:,3,2) = ts_sf(:,:,2,1);
bias_data(:,:,3,3) = tf_sf(:,:,2,1);
bias_data(:,:,3,4) = ti_sf(:,:,2,1);

% feature
bias_data(:,:,4,1) = tt_sf(:,:,2,2);
bias_data(:,:,4,2) = ts_sf(:,:,2,2);
bias_data(:,:,4,3) = tf_sf(:,:,2,2);
bias_data(:,:,4,4) = ti_sf(:,:,2,2);

%% Plot (lots of subpanels)

colorOpts = {[0 0 0],colors.shared,colors.spatial,colors.feature,[0 0 0]};

reportType = {'color','direction'};
params = {'target','side','feature','distractor'};

for cond = 2
    for param = 3
        h = figure; hold on
        
        bdata = squeeze(bias_data(:,cond,:,param));
        
        for pos = 1:5
            plot(pos,bdata(:,pos),'o','MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor','w','MarkerSize',3);
        end
        
        mu = squeeze(nanmean(bdata));
        ci95 = bootci(10000,@nanmean,bdata);
        
        for pos = 1:5
            plot(pos,mu(pos),'o','MarkerFaceColor',colorOpts{pos},'MarkerEdgeColor','w','MarkerSize',5);
        end
        
        axis([0 6 0 1]);
        
        set(gca,'XTick',1:5);
        set(gca,'YTick',[0 1]);
        
        drawPublishAxis('figSize=[2.1,2.5]');
        
        savepdf(h,fullfile('~/proj/afcom/figures/',sprintf('%s_%s.pdf',reportType{cond},params{param})));
    end
end

%% Statistics for bias parameters

models = {'uncued','shared','spatial','feature','no-distractor'};
% check if shared - unshared went up
for cond = 1:2
    for param = 3
        bdata = squeeze(bias_data(:,cond,:,param));
        disp(sprintf('beta%s',params{param}));
        for mi = 2
            disp(sprintf('report-%s %s, model: %s',reportType{cond},params{param},models{mi}));
            disp(sprintf('shared: %1.2f%% %s: %1.2f%%',mean(bdata(:,1))*100,models{mi},100*mean(bdata(:,mi))));
            [p,h,stats] = signrank(bdata(:,mi)-bdata(:,1),0,'tail','left');
            disp(sprintf('P = %1.2f',p));
        end
    end
end

% check if any of the spatial/feature bias parameters are different from the shared
% parametr
% check if bias 3 / 4 are different from 2

for cond = 1:2
    for param = 3:4
        shared = squeeze(bias_data(:,cond,2,:));
        select = squeeze(bias_data(:,cond,param,:));
        % check if any params are different
        disp(sprintf('Checking %s against shared',models{param}));
        for bi = 1:4
            [p,h,stats] = signrank(shared(:,bi),select(:,bi));
            disp(sprintf('%s P=%1.2f',params{bi},p));
        end
    end
end






%% OLD CODE: NOT USED

%% Plot the % of different bias and sensitivity plots



for cond = 1:2
    h = figure;
    for mi = 1:4
        s = subplot(4,1,mi); 

        cbs = bs_(cond,mi);
        cbf = bf_(cond,mi);
        cbi = bi_(cond,mi);
        
        addBiasPlot(s,cbs,cbf,cbi,ac_cmap);
        drawPublishAxis('figSize=[4.4,5]');
        title(models{mi});

    end

    savepdf(h,fullfile('~/proj/afcom/figures/',sprintf('bias_report%s.pdf',reportType{cond})));
end

bs_sf_ = squeeze(mean(bs_sf));
bf_sf_ = squeeze(mean(bf_sf));
bi_sf_ = squeeze(mean(bi_sf));

sf_models = {'shared_sens','shared_all'};

for cond = 1:2
        h = figure;
    for mi = 1:2
        
        select = {'spatial','feature'};
        
        for si = 1:2
            s = subplot(2,2,(mi-1)*2 + si); 

            cbs = bs_sf_(cond,mi,si);
            cbf = bf_sf_(cond,mi,si);
            cbi = bi_sf_(cond,mi,si);
            addBiasPlot(s,cbs,cbf,cbi,ac_cmap);
            drawPublishAxis('figSize=[4.4,5]');
            
            title(sprintf('%s %s',sf_models{mi},select{si}));
        end
    end

    savepdf(h,fullfile('~/proj/afcom/figures/',sprintf('bias_shared_report%s.pdf',reportType{cond})));
end