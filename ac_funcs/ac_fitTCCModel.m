function fit = ac_fitTCCModel(adata,mode,fit)
%% AC_FITTCCMODEL
% This version of the model fits the "TCC" model from this paper:
% Schurgin, M. W., Wixted, J. T., & Brady, T. F. (2018). Psychological scaling reveals a single parameter framework for visual working memory. bioRxiv, 325472.
%
%
%
%
% PARAMETERS
%       guess rate
%   lambda - lapse rate (lambda) and response rate (1-lambda)
%       betas: what % of the time do you pick the other stimuli
%   bs - beta_correctside 
%   bf - beta_correctfeat
%   bi - beta_irrelevant
%       dprimes: how well do you encode the four different stimuli
%   dt - dprime_target
%   ds - dprime_side
%   df - dprime_feature
%   di - dprime_irrelevant
%
% MODE OPTIONS
%
%   you can pass in the conditions that you want to be fit
%   (all,spatial,feature, etc) this allows you to fit the model to
%   different sets of conditions, for example just the "all" condition, or
%   the combination of the spatial/feature conditions, etc. 
%
%   For example the model string:
%     spatial,feature,sh_sens
%   Would fit a model to the cue 2 spatial/feature trials with a single set
%   of a sensitivity parameters shared across conditions
%
%   nosens = per-condition bias parameters, with d' *shared* across conditions
%   nobias = per-condition d' parameters, bias *shared* across conditions
%   nosens,nobias = shared bias and d' parameters
%

global fixedParams

fixedParams.trialTypes = {'all','spatial','feature','target','baseline'};
fixedParams.trialTypeVals = [0 1 2 3 4];
keep = cellfun(@(x) ~isempty(strfind(mode,x)),fixedParams.trialTypes);
fixedParams.trialTypes = fixedParams.trialTypes(keep);
fixedParams.trialTypeVals = fixedParams.trialTypeVals(keep);

% remove trials we aren't fitting
if length(fixedParams.trialTypes)<5
    disp(sprintf('Removing %i conditions that are not being fit',5-length(fixedParams.trialTypes)));
    trialTypes = adata(:,2);
    idxs = zeros(size(trialTypes));
    for tt = 1:length(fixedParams.trialTypeVals)
        idxs = idxs + (trialTypes==fixedParams.trialTypeVals(tt));
    end
    disp(sprintf('Dropped %i trials',size(adata,1)-sum(idxs>0)));
    adata = adata(logical(idxs),:);
end

fixedParams.bdist = ~isempty(strfind(mode,'bdist'));

if ~isempty(strfind(mode,'recovery'))
    fixedParams.trialTypes = {'recovery'};
    fixedParams.trialTypeVals = 0;
end
    
%% TEST FIT
if ~isempty(strfind(mode,'eval'))
    % evaluate fit of existing model
    [~,fit] = vmlike(fit.rawparams,adata,1);
    return
else
    fit = struct;
end

%% CROSS-VALIDATION
if isempty(strfind(mode,'nocv'))
    disp('Running cross-validation');
    % split adata according to runs
    runs = unique(adata(:,15));
    aprobs = [];
    
    for ri = 1:length(runs)
        test = sel(adata,15,runs(ri));
        train = fil(adata,15,'!=',runs(ri));
        
        if ~isempty(intersect(test(:,15),train(:,15)))
            warning('Failure to fully separate training and testing data');
            keyboard
        end
        
        trainfit = ac_fitTCCModel(train,strcat(mode,',nocv'));
        testfit = ac_fitTCCModel(test,strcat(mode,',eval'),trainfit);
        aprobs = [aprobs; testfit.probs];
        like(ri) = testfit.likelihood;
    end
    
    fit = ac_fitTCCModel(adata,strcat(mode,',nocv'));
    fit.cv.probs = aprobs;
    fit.cv.likelihood = -nansum(log(aprobs));
    return
end

%% PARAMETERS

% params.lambda = [0.02 0 1 0 0.1];

% set up any shared parameters
if ~isempty(strfind(mode,'sh_sens'))
    % sensitivity parameter shared for all conditions
    fixedParams.shared_sensitivity = true;
    fixedParams.one_sensitivity = false;
    params.dt_sh = [2 0.01 10 1 3];
else
    fixedParams.shared_sensitivity = false;
    fixedParams.one_sensitivity = false;
    if ~isempty(strfind(mode,'cued_sens'))
        % Separate sensitivity parameters for uncued and cued
        fixedParams.cued_sensitivity = true;
        % we need two sets of parameters, one for uncued and one for cued
        % (shared across both cue types)
        params.dt_1 = [2 0.01 10 1 3];
        % these will be split into dx_2 and dx_3 inside the function
        params.dt_cu = [2 0.01 10 1 3];
    else
        % separate sensitivity parameters for *all* conditions
        fixedParams.cued_sensitivity = false;
        for tt = 1:length(fixedParams.trialTypes)
            params.(sprintf('dt_%i',tt)) = [2 0.01 10 1 3];
        end
    end
end

if ~isempty(strfind(mode,'sh_bias'))
    fixedParams.shared_bias = true;
    params.bs_sh = [0.75 0 1 0.5 1];
    params.bf_sh = [0.75 0 1 0.5 1];
    params.bi_sh = [0.1 0 1 0 0.2];
else
    fixedParams.shared_bias = false;
    if ~isempty(strfind(mode,'cued_bias'))
        fixedParams.cued_bias = true;
        params.bs_1 = [0.75 0 1 0.5 1];
        params.bf_1 = [0.75 0 1 0.5 1];
        params.bi_1 = [0.75 0 1 0.5 1];
        
        params.bs_cu = [0.75 0 1 0.5 1];
        params.bf_cu = [0.75 0 1 0.5 1];
        params.bi_cu = [0.75 0 1 0.5 1];
    else
        fixedParams.cued_bias = false;
        for tt = 1:length(fixedParams.trialTypes)
            params.(sprintf('bs_%i',tt)) = [0.75 0 1 0.5 1];
            params.(sprintf('bf_%i',tt)) = [0.75 0 1 0.5 1];
            params.(sprintf('bi_%i',tt)) = [0.1 0 1 0 0.2];
        end
    end
end

% don't fit a bias unless requested (not useful, nobody has a bias)
% if ~isempty(strfind(mode,'bias'))
%     params.bias = [0 -pi pi -pi/8 pi/8];
% else
%     params.bias = 0; % [0 -pi pi];
% end

[ip,minp,maxp,plb,pub] = initParams(params);

% BADS VERSION
if strfind(mode,'bads')
%     warning('Tolerance size is large and maxfunevals is tiny: reduce for main fits');
    options.TolMesh = 0.001;
%     options.MaxFunEvals = 200;
    
    bestparams = bads(@(p) vmlike(p,adata,0),ip,minp,maxp,plb,pub,[],options);

% FMINCON VERSION
else
    warning('Tolerance size is large: reduce for main fits');
%     options = optimoptions('fmincon','Algorithm','active-set','TolFun',5,'TolCon',1,'Display','off'); % set a limit or it goes on foreeeeeeeeeeeever

    bestparams = fmincon(@(p) vmlike(p,adata,0),ip,[],[],[],[],minp,maxp,[]);%,options);
end

[~,fit] = vmlike(bestparams,adata,1);

fit.rawparams = bestparams;
fit.trialTypes = fixedParams.trialTypes;
fit.trialTypeVals = fixedParams.trialTypeVals;
fit.data = adata;

function [likelihood, fit] = vmlike(params,adata,computeOutput)
global fixedParams


if any(isnan(params))
    warning('parameter evaluated to nan');
    stop = 1;
    likelihood = inf; 
    return
end

params = getParams(params);

% headers = 
% 
%   Columns 1 through 6
% 
%     'target'    'trialType'    'cue'    'duration'    'dead'    'targetAngle'
% 
%   Columns 7 through 12
% 
%     'distractorAngle'    'angle1'    'angle2'    'angle3'    'angle4'    'respAngle'
% 
%   Columns 13 through 14
% 
%     'respDistance'    'distDistance'
% 


% setup probabilities list
probs = zeros(size(adata,1),1);

% pdft = computePDF(params.

%% 
% test code
% xs = -pi:pi/64:pi;
% dprimes = 0.2:0.5:3;
% figure; hold on
% legs = {};
% clear out outs
% for di = 1:length(dprimes)
%     dprime = dprimes(di);
%     tic
%     out = preComputeTCCPDF(xs,dprime);
%     toc
%     outs{di} = out;
%     sum(outs{di})
%     plot(xs,out);
%     legs{end+1} = sprintf('%1.2f',dprime);
% end
% legend(legs);

%% Compute for each dprime parameter the likelihood functions
xs = -pi:pi/256:pi;
for tt = 1:length(fixedParams.trialTypes)
    
    idxs = adata(:,2)==fixedParams.trialTypeVals(tt);
    % get trials with this trial type
    tdata = adata(idxs,:);
    % for each trial, re-organize the angles so we can get the likelihoods
    % each column is the target, side, feat, and irrelevant angles all in
    % one calculation
    angleOpts = [8     9    10    11
                9     8    11    10
                10    11     8     9
                11    10     9     8];

    angles = zeros(size(tdata,1),4);
    for ti = 1:size(tdata,1)
        trial = tdata(ti,:);
        angles(ti,:) = trial(angleOpts(trial(1),:));
    end

    % rotate all the angles relative to the response angle 
    rt = repmat(tdata(:,12),1,4);
    angles = sign(angles-rt).*angdist(rt,angles);
    
    % is that right? or should it be:
%     angles = mod(angles-repmat(tdata(:,12),1,4),2*pi);

    % get the parameters for this trialtype
    if fixedParams.one_sensitivity && fixedParams.shared_sensitivity
        dt = params.d_sh;
        ds = params.d_sh;
        df = params.d_sh;
        di = params.d_sh;
    elseif fixedParams.shared_sensitivity
        dt = params.dt_sh;
        ds = params.dt_sh;
        df = params.dt_sh;
        di = params.dt_sh;
    else
        if fixedParams.cued_sensitivity
            if tt==1
                dt = params.dt_1;
                ds = params.dt_1;
                df = params.dt_1;
                di = params.dt_1;
            else
                dt = params.dt_cu;
                ds = params.dt_cu;
                df = params.dt_cu;
                di = params.dt_cu;
            end
        else
            dt = params.(sprintf('dt_%i',tt));
            ds = params.(sprintf('dt_%i',tt));
            df = params.(sprintf('dt_%i',tt));
            di = params.(sprintf('dt_%i',tt));
        end
    end
    
    liket = preComputeTCCPDF(xs,dt);
%     likes = preComputeTCCPDF(xs,ds);
%     likef = preComputeTCCPDF(xs,df);
%     likei = preComputeTCCPDF(xs,di);

    like = zeros(size(angles));
    % compute the probability that the response angle was pulled from this
    % distribution
    like(:,1) = interp1(xs,liket,angles(:,1),'linear');
    like(:,2) = interp1(xs,liket,angles(:,2),'linear');
    like(:,3) = interp1(xs,liket,angles(:,3),'linear');
    like(:,4) = interp1(xs,liket,angles(:,4),'linear');

    % weight the likelihoods by the beta values
    if fixedParams.shared_bias
        bs = params.bs_sh;
        bf = params.bf_sh;
        bi = params.bi_sh;
    else
        if fixedParams.cued_bias
            if tt==1
                bs = params.bs_1;
                bf = params.bf_1;
                bi = params.bi_1;
            else
                bs = params.bs_cu;
                bf = params.bf_cu;
                bi = params.bi_cu;
            end
        else
            bs = params.(sprintf('bs_%i',tt));
            bf = params.(sprintf('bf_%i',tt));
            bi = params.(sprintf('bi_%i',tt));
        end
    end
    betas = [bs*bf bs*(1-bf) (1-bs)*(1-bi) (1-bs)*bi];

    trial_likelihoods = like * betas';
    
    probs(idxs) = trial_likelihoods;

    if any(trial_likelihoods==0)
        stop = 1;
    end
end

likelihood = -nansum(log(probs));

if isinf(likelihood)
    warning('Likelihood evaluated to Inf');
    probs(probs==0) = eps;
    likelihood = -nansum(log(probs));
end

if isnan(likelihood)
    warning('Likelihood evaluated to NaN');
    stop = 1;
end

fit.probs = probs;
fit.likelihood = likelihood;
fit.params = params;

if computeOutput
    
    x = -pi:pi/128:pi;
    fit.x = x;
    fit.trialTypes = fixedParams.trialTypes;
    fit.targetOrder = {'target','side','feature','distractor'};
    
    out = zeros(5,4,length(x));
    outs = out;
    % for each trial type, compute the likelihood function?
    tx = 0:pi/128:pi;
    for tt = 1:length(fixedParams.trialTypes)
        
        if fixedParams.one_sensitivity && fixedParams.shared_sensitivity
            dt = params.d_sh;
            ds = params.d_sh;
            df = params.d_sh;
            di = params.d_sh;
        elseif fixedParams.shared_sensitivity
            dt = params.dt_sh;
            ds = params.dt_sh;
            df = params.dt_sh;
            di = params.dt_sh;
        else
            if fixedParams.cued_sensitivity
                if tt==1
                    dt = params.dt_1;
                    ds = params.dt_1;
                    df = params.dt_1;
                    di = params.dt_1;
                else
                    dt = params.dt_cu;
                    ds = params.dt_cu;
                    df = params.dt_cu;
                    di = params.dt_cu;
                end
            else
                dt = params.(sprintf('dt_%i',tt));
                ds = params.(sprintf('dt_%i',tt));
                df = params.(sprintf('dt_%i',tt));
                di = params.(sprintf('dt_%i',tt));
            end
        end

        liket = computeTCCPDF(tx,dt);
        liket = [fliplr(liket) liket(2:end)];
        likes = computeTCCPDF(tx,ds);
        likes = [fliplr(likes) likes(2:end)];
        likef = computeTCCPDF(tx,df);
        likef = [fliplr(likef) likef(2:end)];
        likei = computeTCCPDF(tx,di);
        likei = [fliplr(likei) likei(2:end)];

        % normalize everything
        liket = liket./sum(liket);
        likes = likes./sum(likes);
        likef = likef./sum(likef);
        likei = likei./sum(likei);

        if fixedParams.shared_bias
            bs = params.bs_sh;
            bf = params.bf_sh;
            bi = params.bi_sh;
        else
            if fixedParams.cued_bias
                if tt==1
                    bs = params.bs_1;
                    bf = params.bf_1;
                    bi = params.bi_1;
                else
                    bs = params.bs_cu;
                    bf = params.bf_cu;
                    bi = params.bi_cu;
                end
            else
                bs = params.(sprintf('bs_%i',tt));
                bf = params.(sprintf('bf_%i',tt));
                bi = params.(sprintf('bi_%i',tt));
            end
        end
        % scale the likelihood functions 

        % target/side/feat/dist
        out(tt,1,:) = liket;
        out(tt,2,:) = likes;
        out(tt,3,:) = likef;
        out(tt,4,:) = likei;

        outs(tt,1,:) = bs*bf*liket;
        outs(tt,2,:) = bs*(1-bf)*likes;
        outs(tt,3,:) = (1-bs)*(1-bi)*likef;
        outs(tt,4,:) = (1-bs)*bi*likei;
    end    
    fit.outs = outs;
    fit.out = out;
end

%% Helper routines

function d = angdist(t1,t2)
d = acos(cos(t1).*cos(t2)+sin(t1).*sin(t2));
%%
function [initparams, minparams, maxparams, plb, pub] = initParams(params)

global fixedParams

fixedParams.strs = fields(params);
fixedParams.num = length(fixedParams.strs);

initparams = [];
minparams = [];
maxparams = [];
plb = [];
pub = [];
indexes = zeros(1,fixedParams.num);
count = 1;

fixed = zeros(1,fixedParams.num);
optim = zeros(1,fixedParams.num);

for i = 1:fixedParams.num
    cvals = params.(fixedParams.strs{i});
    
    if length(cvals)==1
        fixedParams.(fixedParams.strs{i}) = cvals;
        fixed(i) = 1;
    elseif length(cvals)==3
        initparams = [initparams cvals(1)];
        minparams = [minparams cvals(2)];
        maxparams = [maxparams cvals(3)];
        plb = [plb cvals(2)];
        pub = [pub cvals(3)];
        indexes(i) = count;
        count = count+1;
    elseif length(cvals)==5
        initparams = [initparams cvals(1)];
        minparams = [minparams cvals(2)];
        maxparams = [maxparams cvals(3)];
        plb = [plb cvals(4)];
        pub = [pub cvals(5)];
        indexes(i) = count;
        count = count+1;
    elseif length(cvals)==2 || length(cvals)>3
        % optimizer
        fixedParams.(fixedParams.strs{i}) = cvals;
        optim(i) = 1;
    else
        error('You initialized a parameter with the wrong initial values... unable to interpret');
    end
end
fixedParams.optim = optim;
fixedParams.fixed = fixed;
fixedParams.idx = indexes;

function p = getParams(params)

global fixedParams

for i = 1:fixedParams.num
    if fixedParams.fixed(i)
        p.(fixedParams.strs{i}) = fixedParams.(fixedParams.strs{i});
    else
        p.(fixedParams.strs{i}) = params(fixedParams.idx(i));
    end
end