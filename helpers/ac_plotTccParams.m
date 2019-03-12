function ac_plotTccParams(fit)

%%

cmap = ac_cmap;

params = fit.params;
trialTypes = fit.trialTypes;
trialTypeVals = fit.trialTypeVals;

shared_axis = [-180 180 0 0.015];

groups = {'cue 4','cue 2 spatial','cue 2 feature','cue 1 target','cue 1 baseline'};
groups = {'cue_4','cue_side','cue_feat','cue_1','cue_target'};

x = -180:180;
p = 1-pscale(abs(x));

for tt = 1:length(trialTypes)
    h = figure;
    % pull the betas and dprimes
    bs = params.(sprintf('bs_%i',tt));
    bf = params.(sprintf('bf_%i',tt));
    bi = params.(sprintf('bi_%i',tt));
    
    dt = params.(sprintf('dt_%i',tt));
    df = params.(sprintf('df_%i',tt));
    ds = params.(sprintf('ds_%i',tt));
    di = params.(sprintf('di_%i',tt));
    
    % come up with a way to display this? 
    
    % first compute the dprime curves
%     for xi = 1:length(x)
%         pt(xi) = sum(log(normcdf((dt*(p(xi)-p)))));
%         pf(xi) = sum(log(normcdf((df*(p(xi)-p)))));
%         ps(xi) = sum(log(normcdf((ds*(p(xi)-p)))));
%         pi(xi) = sum(log(normcdf((di*(p(xi)-p)))));
%     end

    pt = dt*p;
    pf = df*p;
    ps = ds*p;
    pi = di*p;
    
    pt = pt./sum(pt);
    pf = pf./sum(pf);
    ps = ps./sum(ps);
    pi = pi./sum(pi);
    
    % weight by the right scaling
    subplot(141);
    plot(x,pt * bs * bf,'Color',cmap.target);
    axis(shared_axis);
    subplot(142);
    plot(x,pf * bs * (1-bf),'Color',cmap.side);
    axis(shared_axis);
    subplot(143);
    plot(x,ps * (1-bs) * (1-bf),'Color',cmap.feat);
    axis(shared_axis);
    subplot(144);
    plot(x,pi * (1-bs) * bi,'Color',cmap.dist);
    axis(shared_axis);
    
    savepdf(h,fullfile('~/proj/afcom/figures',fit.dataType,sprintf('TCC_%s',groups{tt})));
end