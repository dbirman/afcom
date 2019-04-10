function ac_plotTccParams(fit)

%%

cmap = ac_cmap;

params = fit.params;
trialTypes = fit.trialTypes;
trialTypeVals = fit.trialTypeVals;

shared_axis = [-pi pi 0 0.015];

groups = {'cue 4','cue 2 spatial','cue 2 feature','cue 1 target','cue 1 baseline'};
groups = {'cue_4','cue_side','cue_feat','cue_1','cue_target'};

x = -pi:.01:pi;

% first compute the baseline
tt = 5;
[lt_t,lf_t,ls_t,li_t] = getLikeFromParams(x,params,tt);

for tt = 1:(length(trialTypes)-1)
    h = figure;
    
    [lt,lf,ls,li] = getLikeFromParams(x,params,tt);
    
    % weight by the right scaling
    subplot(141); hold on
    plot(x,lt_t,'Color','k');
    plot(x,lt,'Color',cmap.target);
    axis(shared_axis);
    subplot(142); hold on
    plot(x,lf_t,'Color','k');
    plot(x,lf,'Color',cmap.side);
    axis(shared_axis);
    subplot(143); hold on
    plot(x,ls_t,'Color','k');
    plot(x,ls,'Color',cmap.feat);
    axis(shared_axis);
    subplot(144); hold on
    plot(x,li_t,'Color','k');
    plot(x,li,'Color',cmap.dist);
    axis(shared_axis);
    
    savepdf(h,fullfile('~/proj/afcom/figures',fit.dataType,sprintf('TCC_%s',groups{tt})));
end

function [lt,lf,ls,li] = getLikeFromParams(x,params,tt)

% pull the betas and dprimes
bs = params.(sprintf('bs_%i',tt));
bf = params.(sprintf('bf_%i',tt));
bi = params.(sprintf('bi_%i',tt));

dt = params.(sprintf('dt_%i',tt));
df = params.(sprintf('df_%i',tt));
ds = params.(sprintf('ds_%i',tt));
di = params.(sprintf('di_%i',tt));

% come up with a way to display this? 

lt = computeTCCPDF(x,dt);
lt = lt*bs*bf;
lf = computeTCCPDF(x,df);
lf = lf*bs*(1-bf);
ls = computeTCCPDF(x,ds);
ls = ls*(1-bs)*(1-bi);
li = computeTCCPDF(x,di);
li = li*(1-bs)*bi;