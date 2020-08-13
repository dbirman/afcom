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

for dpinc = [repmat(0.1,1,100) repmat(0.25,1,100) repmat(0.5,1,100) repmat(1,1,100)]
    cuedp = basedp+dpinc;
    cueLike = computeTCCPDF(-pi:pi/128:pi,cuedp);
    
    
    adata = nan(1010,16);
    
    for n = 1:10
        % uncued trials (tt==1)
        offAngles = rand(1,3)*2*pi;
        
        % generate the likelihood functions
        offLike = nan(length(x),4);
        offLike(:,1) = baseLike;
        for oi = 1:3
            shiftIdx = find((x+pi)>=offAngles(oi),1);
            offLike(:,oi+1) = circshift(baseLike,[0,shiftIdx]);
        end
        
        finalLike = cumsum(offLike*betas');
        finalLike = finalLike./max(finalLike); % normalize in case the likelihood is messed up
        
        resp = x(find(finalLike>=rand,1))+pi;
        
        adata(n,:) = [1 0 nan 1 0 pi nan pi offAngles resp nan nan nan nan];
    end
    
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
        
        adata(n+500,:) = [1 1 nan 1 0 pi nan pi offAngles resp nan nan nan nan];
    end
    
    dataset{end+1} = adata;
end

%% Fit models
clear fit
parfor di = 1:length(dataset)
    fit{di} = ac_fitTCCModel(dataset{di},'bads,all,spatial,feature,cued_sens,sh_bias,nocv');
end

save(fullfile('~/proj/afcom/dprime_recovery.mat'),'dataset','fit');

%% Figures
load(fullfile('~/proj/afcom/dprime_recovery.mat'));

exp = [repmat(0.1,1,100) repmat(0.25,1,100) repmat(0.5,1,100) repmat(1,1,100)];
for fi = 1:400
    p = fit{fi}.params;
    
    db(fi) = p.dt_1;
    dc(fi) = p.dt_cu;
end

%% Figure 1
exp = 1+[0.1 0.25 0.5 1];
clear mu ci
for i = 1:4
    vals = dc((((i-1)*100)+1):((i-1)*100+100));
    mu(i) = mean(vals);
    ci(i) = 1.96*std(vals);
end

h = figure; hold on
errbar(exp,mu,ci,'-k');
plot(exp,mu,'o','MarkerFaceColor','k','MarkerEdgeColor','w','MarkerSize',5);
axis([1 2 1 2]);
plot([1 2],[1 2],'--k');
set(gca,'XTick',[1 2],'YTick',[1 2]);
xlabel('Simulated d');
ylabel('Fitted d');

drawPublishAxis('figSize=[4.5, 2.5]');

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
    mu(g,:) = mean(cbetas);
    ci(g,:) = std(cbetas);
%     ci(g,:,:) = bootci(1000,@nanmean,cbetas);
end

%% Display the probability of recovering different values 
% rng = squeeze(ci(:,2,:))-mu;
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
    drawPublishAxis('figSize=[4.5,6.5]');
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