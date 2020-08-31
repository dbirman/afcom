%%

% TARGET TRIALTYPE CUE DURATION DEAD TARGETANGLE DISTRACTORANGLE
%    1      2       3     4       5       6           7
% ANGLE1 ANGLE2 ANGLE3 ANGLE4 RESPANGLE RESPDISTANCE DISTDISTANCE RUNS RT
%    8      9     10     11      12          13            14      15  16
   
% Generate fake data with various d' or bias, then try to recover the
% parameters. We'll generate 1000 trials.

% We will test first whether two d' differences can be detected

%% Pure d' test
basedp = 1;

dataset = {};

x = -pi:pi/128:pi;
baseLike = computeTCCPDF(x,basedp);

betas = [1 0 0 0];

odps = [0 0.01 .02 .03 .04 0.05 .1 .2];
reps = 200;

dps = repmat(odps,reps,1);
dps = dps(:)';

for dpinc = dps
    cuedp = basedp+dpinc;
    cueLike = computeTCCPDF(-pi:pi/128:pi,cuedp);
    
    
    adata = nan(700,16);
    
    for n = 1:700
        % cued trials (tt==2)
        offAngles = rand(1,3)*2*pi;
        
        % generate the likelihood functions
        offLike = nan(length(x),4);
        offLike(:,1) = cueLike;
        for oi = 1:3
            shiftIdx = find((x+pi)>=offAngles(oi),1);
            offLike(:,oi+1) = circshift(cueLike,[0,shiftIdx]);
        end
        
        finalLike = cumsum(offLike*betas');
        finalLike = finalLike./max(finalLike); % normalize in case the likelihood is messed up
        
        resp = x(find(finalLike>=rand,1))+pi;
        
        adata(n,:) = [1 0 nan 1 0 pi nan pi offAngles resp nan nan nan nan];
    end
    
    dataset{end+1} = adata;
end

%% Fit models
clear fit
parfor di = 1:length(dataset)
    fit{di} = ac_fitTCCModel(dataset{di},'bads,all,sh_sens,sh_bias,nocv');
end

save(fullfile('~/proj/afcom/dprime_recovery.mat'),'dataset','fit');

%% Now we're going to compute power on these
% i.e. true positive: what percent did we get a result
% vs. false positive: what percent would have been result anyway
load(fullfile('~/proj/afcom/dprime_recovery.mat'));

fitd = zeros(size(dps));
for di = 1:length(dps)
    fitd(di) = fit{di}.params.dt_sh;
end

%% Do the BOOT
dps_ = zeros(8,200);
for di = 1:8
    idx = (di-1)*200;
    dps_(di,:) = fitd(idx+1:idx+200);
end

reps = 10000;
higher = zeros(7,reps);
for di = 2:8
    uncued = randsample(dps_(1,:),reps,true);
    cued = randsample(dps_(di,:),reps,true);
    higher(di-1,:) = cued>uncued;
end

mean(higher,2);

%% Figure
h = figure(1); clf
hold on

higher_ = [0.5 mean(higher,2)'];


% fit model
fit = ac_fitSatExp([odps 0.35],[higher_ 1]);

plot(fit.x_,fit.y_,'-k');
plot(odps,higher_,'o','MarkerFaceColor','k','MarkerEdgeColor','w');
vline(0.11,'-r');
vline(0.33,'-r');
xlabel('d'' effect size');
ylabel('Sensitivity');
drawPublishAxis('figSize=[4,4]');

savepdf(h,fullfile('~/proj/afcom/figures/model_recovery_d.pdf'));

%% BETA TEST
clear all uncued cued

uncued_orig = [0.5 0.35 0.1 0.05];
cued_orig = [0.75 0.25 0 0];
% approximate values

basedp = 1;

reps = 200;

percs = 0:.1:1;
dataset = cell(length(percs),reps);

x = -pi:pi/128:pi;
baseLike = computeTCCPDF(x,basedp);

% let's compute the sensitivity for interpolations between these 
diff = cued-uncued;

for pii = 1:length(percs)
    perc = percs(pii);
    cbeta = uncued + diff*perc;
    disp(cbeta);
    disp(sum(cbeta));
    
    for rep = 1:reps
        adata = nan(700,16);
        for n = 1:700
            % uncued trials (tt==1)
            offAngles = rand(1,3)*2*pi;

            % generate the likelihood functions
            offLike = nan(length(x),4);
            offLike(:,1) = baseLike;
            for oi = 1:3
                shiftIdx = find((x+pi)>=offAngles(oi),1);
                offLike(:,oi+1) = circshift(baseLike,[0,shiftIdx]);
            end

            % do this the way it's supposed to work, use beta to choose which
            % thing we will run off of and then pull the location from there
            choice = find(rand<=cumsum(cbeta),1);
            finalLike = cumsum(offLike(:,choice));
            finalLike = finalLike ./ max(finalLike);

    %         finalLike = cumsum(offLike*cbeta');
    %         finalLike = finalLike./max(finalLike); % normalize in case the likelihood is messed up

            resp = x(find(finalLike>=rand,1));
            if choice==1
                resp = resp+pi;
            end
    %         disp(resp)

            adata(n,:) = [1 0 nan 1 0 pi nan pi offAngles resp nan nan nan nan];
        end
        dataset{pii,rep} = adata;
    end
    
    disp(pii);
end

%% Fit

% best to fit everything
for pii = 1:size(dataset,1)
    for ri = 1:200
        sz((pii-1)*200+ri) = size(dataset{pii,ri},1);
        dataset_{(pii-1)*200+ri} = dataset{pii,ri};
    end
end

clear fit_ fit
parfor di = 1:length(dataset_)
    fit_{di} = ac_fitTCCModel(dataset_{di},'bads,all,spatial,feature,sh_sens,sh_bias,nocv');
end

% reorganize
for pii = 1:size(dataset,1)
    for ri = 1:200
        fit{pii,ri} = fit_{(pii-1)*200+ri};
    end
end

save(fullfile('~/proj/afcom/beta_uncued_cued_recovery.mat'),'dataset','fit');

%% Re-organize into the 4 main betas

load(fullfile('~/proj/afcom/beta_uncued_cued_recovery.mat'));

betas = nan(11,200,4);

for pi = 1:11
    for ri = 1:200
        p = fit{pi,ri}.params;
        
        % compute the betas
        betas(pi,ri,:) = [p.bs_sh*p.bf_sh p.bs_sh*(1-p.bf_sh) (1-p.bs_sh)*(1-p.bi_sh) (1-p.bs_sh)*p.bi_sh];

    end
end

%% Do the BOOT
for bi = 1:4
    clear higher
    beta = betas(:,:,bi);
    dir = [1 -1 -1 -1];

    reps = 10000;
    % we'll ask four questions: did beta(1) go up, and did beta(2-4) go down
    higher = zeros(10,reps);

    % percs(3) is 0, so we'll start from 3 and go up
    for i = 2:11
        uncued = beta(1,randsample(1:200,reps,true));
        cued = beta(i,randsample(1:200,reps,true));
        if bi==1
            higher(i-1,:) = cued > uncued;
        else
            higher(i-1,:) = cued < uncued;
        end
    end

    mean(higher,2)
    
    %%
    h = figure(1); clf
    hold on

    higher_ = [0.5 mean(higher,2)'];


    % fit model
    xs = cued_orig(bi) + percs*diff(bi);
    fit = ac_fitSatExpBeta([percs],[higher_]);

    plot(fit.x_,fit.y_,'-k');
    plot(percs,higher_,'o','MarkerFaceColor','k','MarkerEdgeColor','w');
    
    axis([0 1 0.5 1]);
    set(gca,'XTick',0:.25:1,'YTick',[0.5:.25:1]);
    vline(1,'-r');
% %     vline(0.33,'-r');
    xlabel('Percentage of reported effect');
    ylabel('Sensitivity');
    drawPublishAxis('figSize=[4,4]');

    savepdf(h,fullfile('~/proj/afcom/figures/',sprintf('model_recovery_beta%i.pdf',bi)));
end

%% Figure
h = figure(1); clf
hold on

higher_ = [0.5 mean(higher,2)'];


% fit model
fit = ac_fitSatExp([odps 0.35],[higher_ 1]);

plot(fit.x_,fit.y_,'-k');
plot(odps,higher_,'o','MarkerFaceColor','k','MarkerEdgeColor','w');
vline(0.11,'-r');
vline(0.33,'-r');
xlabel('d'' effect size');
ylabel('Sensitivity');
drawPublishAxis('figSize=[4,4]');

savepdf(h,fullfile('~/proj/afcom/figures/model_recovery_d.pdf'));
%% Plot expected vs actual

%% Pure beta test
% That first test worked! So now let's adjust the beta values without
% adjusting d'
basedp = 1;

dataset = {};

x = -pi:pi/128:pi;
baseLike = computeTCCPDF(x,basedp);

betas_ = [1 0 0 0
         0.5 0.5 0 0
         0.25 0.25 0.25 0.25
         0 0 0 1
         0 0 1 0
         0 1 0 0];
     
betas = [];
for bi = 1:size(betas_,1)
    for rep = 1:100
        betas((bi-1)*100+rep,:) = betas_(bi,:);
    end
end

for bi = 1:size(betas,1)
    cbeta = betas(bi,:);    
    
    adata = nan(1000,16);
    
    for n = 1:700
        % uncued trials (tt==1)
        offAngles = rand(1,3)*2*pi;
        
        % generate the likelihood functions
        offLike = nan(length(x),4);
        offLike(:,1) = baseLike;
        for oi = 1:3
            shiftIdx = find((x+pi)>=offAngles(oi),1);
            offLike(:,oi+1) = circshift(baseLike,[0,shiftIdx]);
        end
        
        % do this the way it's supposed to work, use beta to choose which
        % thing we will run off of and then pull the location from there
        choice = find(rand<=cumsum(cbeta),1);
        finalLike = cumsum(offLike(:,choice));
        finalLike = finalLike ./ max(finalLike);
        
%         finalLike = cumsum(offLike*cbeta');
%         finalLike = finalLike./max(finalLike); % normalize in case the likelihood is messed up
        
        resp = x(find(finalLike>=rand,1));
        if choice==1
            resp = resp+pi;
        end
%         disp(resp)
        
        adata(n,:) = [1 0 nan 1 0 pi nan pi offAngles resp nan nan nan nan];
    end
    
    dataset{end+1} = adata;
end

%% Fit models
clear fit
parfor di = 1:length(dataset)
    fit{di} = ac_fitTCCModel(dataset{di},'bads,all,spatial,feature,sh_sens,sh_bias,nocv');
end

save(fullfile('~/proj/afcom/beta_recovery.mat'),'dataset','fit');
%%
for di = 1:length(fit)
    p = fit{di}.params;
    betas = [p.bs_sh*p.bf_sh p.bs_sh*(1-p.bf_sh) (1-p.bs_sh)*(1-p.bi_sh) (1-p.bs_sh)*p.bi_sh];
    disp(sprintf('%1.2f %1.2f %1.2f %1.2f',round(betas*100)/100));
%     disp(round(betas))
end

%% 
load(fullfile('~/proj/afcom/beta_recovery.mat'));
clear betas
for i = 1:600
    p = fit{i}.params;
    betas(i,:) = [p.bs_sh*p.bf_sh p.bs_sh*(1-p.bf_sh) (1-p.bs_sh)*(1-p.bi_sh) (1-p.bs_sh)*p.bi_sh];
end

clear mu ci
for g = 1:6
    cbetas = betas(((g-1)*100+1):g*100,:);
    clear mu_
    for bi = 1:4
        for r = 1:1000
            sample = randsample(cbetas(:,bi),7);
            mu_(r) = mean(sample);
        end
        mu(g,bi) = mean(mu_);
        ci(g,bi) = std(mu_);
    end
end

%% Display the probability of recovering different values 
rng = ci;
h = figure; hold on
plot([0 1],[0 1],'--k');
errbar(betas_(:),mu(:),rng(:),'-k');
plot(betas_(:),mu(:),'o','MarkerFaceColor','k','MarkerEdgeColor','w','MarkerSize',5);

set(gca,'XTick',[0 1],'YTick',[0 1]);
xlabel('Simulated parameter value');
ylabel('Fitted value');
drawPublishAxis('figSize=[4.5,2.5]');
savepdf(h,fullfile('~/proj/afcom/figures/model_recovery_beta.pdf'));
%% Figure

h = figure(1); clf
for i = 1:6
    subplot(6,1,i); hold on
    
    plot(1:4,betas_(i,:),'--r');
    errbar(1:4,mu(i,:),rng(i,:),'-','Color','k');
    plot(1:4,mu(i,:),'o','MarkerFaceColor','k','MarkerEdgeColor','w','MarkerSize',5);
    
    set(gca,'XTick',1:4,'XTickLabel',{'Target','Side','Feature','Distractor'});
    set(gca,'YTick',[0 1]);
    axis([1 4 0 1]);
    drawPublishAxis('figSize=[4.5,8.9]');
end

savepdf(h,fullfile('~/proj/afcom/figures/model_recovery_varbeta.pdf'));

%% Do some statistics


%% Weird drop in d' tests

clear fit 

basedp = 1;

dataset = {};

x = -pi:pi/128:pi;
baseLike = computeTCCPDF(x,basedp);

betas = [1 0 0 0];

clear adatas
for rep = 1:100
    adata = nan(10,16);

    for n = 1:10
        % uncued trials (tt==1)

        finalLike = cumsum(baseLike);
        finalLike = finalLike./max(finalLike); % normalize in case the likelihood is messed up

        resp = x(find(finalLike>=rand,1))+pi;

        adata(n,:) = [1 0 nan 1 0 pi nan pi offAngles resp nan nan nan nan];
    end
    
    adatas{rep} = adata;
end

parfor rep = 1:100
    fit{rep} = ac_fitTCCModel(adatas{rep},'bads,all,spatial,feature,cued_sens,sh_bias,nocv');
end

%% Figure
clear d
for rep = 1:100
    p = fit{rep}.params;
    d(rep) = p.dt_1;
end

%% Okay that stuff all worked
% Now we're going to set up model recovery for the *actual* results in the
% paper:
% (1) Difference in d' with all betas the same
% (2) Difference in betas with d' the same
% (3) Difference in d' and beta (uncued vs. cued)
% (4) Difference in d' and beta (uncued vs. spatial vs. feature)

