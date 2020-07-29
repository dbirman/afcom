function fit = ac_fitTCCModel_corrected(adata,mode,fit)
%% AC_FITTCCMODEL
% This version of the model fits the "TCC" model from this paper:
% Schurgin, M. W., Wixted, J. T., & Brady, T. F. (2018). Psychological scaling reveals a single parameter framework for visual working memory. bioRxiv, 325472.
%
% In this variation (corrected) the model fits the full likelihood function
% for *each* trial. That means computing the integral each time the model
% fit is run. Unfortunately this makes the model fitting procedure very
% slow compared to the faster version implemented in ac_fitTCCModel
% originally. 
%
% In this variation the beta parameters encode the relative the amplitudes
% against the target d' for how high the likelihoods are. They encode a
% percentage of the d' in this way. So the channel response for c_theta
% ends up being 1*c_theta,target + bs*c_theta,side + bf*c_theta,feat +
% bi*c_theta_irrelevant. Because this happens before the maximum value
% computation the model is "different" from the version where this is all
% computed as sampling from the distributions. Whether the model fit is
% different or not is a matter of comparison. We'll see how slow it is!
%
% PARAMETERS
%       betas: what % of the time do you pick the other stimuli
%   bs - beta_correctside 
%   bf - beta_correctfeat
%   bi - beta_irrelevant
%       dprime: how well do you encode the four different stimuli
%   d
%
% MODE OPTIONS
%
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
    params.dt_sh = [1.5 0.01 10 0.75 3];
else
    fixedParams.shared_sensitivity = false;
    fixedParams.one_sensitivity = false;
    if ~isempty(strfind(mode,'cued_sens'))
        % Separate sensitivity parameters for uncued and cued
        fixedParams.cued_sensitivity = true;
        % we need two sets of parameters, one for uncued and one for cued
        % (shared across both cue types)
        params.dt_1 = [1.5 0.01 10 0.75 3];
        % these will be split into dx_2 and dx_3 inside the function
        params.dt_cu = [1.5 0.01 10 0.75 3];
    else
        % separate sensitivity parameters for *all* conditions
        fixedParams.cued_sensitivity = false;
        for tt = 1:length(fixedParams.trialTypes)
            params.(sprintf('dt_%i',tt)) = [1.5 0.01 10 0.75 3];
        end
    end
end

if ~isempty(strfind(mode,'sh_bias'))
    fixedParams.shared_bias = true;
    params.bs_sh = [1 0.01 10 0.25 3];
    params.bf_sh = [1 0.01 10 0.25 3];
    params.bi_sh = [0.5 0 10 0.25 3];
else
    fixedParams.shared_bias = false;
    if ~isempty(strfind(mode,'cued_bias'))
        fixedParams.cued_bias = true;
        params.bs_1 = [1 0.01 10 0.25 3];
        params.bf_1 = [1 0.01 10 0.25 3];
        params.bi_1 = [0.5 0.01 10 0.25 3];
        
        params.bs_cu = [1 0.01 10 0.25 3];
        params.bf_cu = [1 0.01 10 0.25 3];
        params.bi_cu = [0.5 0.01 10 0.25 3];
    else
        fixedParams.cued_bias = false;
        for tt = 1:length(fixedParams.trialTypes)
            params.(sprintf('bs_%i',tt)) = [0.5 0 10 0 3];
            params.(sprintf('bf_%i',tt)) = [0.5 0 10 0 3];
            params.(sprintf('bi_%i',tt)) = [0.25 0 10 0 3];
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
    warning('Tolerance size is large and maxfunevals is tiny: reduce for main fits');
%     options.TolMesh = 0.001;
    options.TolMesh = 1;
    options.MaxFunEvals = 50;
    
    bestparams = bads(@(p) vmlike(p,adata,0),ip,minp,maxp,plb,pub,[],options);

% FMINCON VERSION
else
    warning('Tolerance size is large: reduce for main fits');
    options = optimoptions('fmincon','Algorithm','active-set','TolFun',5,'TolCon',1,'Display','off'); % set a limit or it goes on foreeeeeeeeeeeever

    bestparams = fmincon(@(p) vmlike(p,adata,0),ip,[],[],[],[],minp,maxp,[],options);
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

% get the predicted likelihood for each response relative to the target
% angle
probs = nan(size(adata,1),1);

%% Build the channel model
% For each trial, get the channel responses
N = 100; % number of trials

% set up the encoder response-range


% for each trial, re-organize the angles so we can get the likelihoods
% each column is the target, side, feat, and irrelevant angles all in
% one calculation
angleOpts = [8     9    10    11
                9     8    11    10
                10    11     8     9
                11    10     9     8];

angles = zeros(size(adata,1),4);
for ti = 1:size(adata,1)
    trial = adata(ti,:);
    target = trial(1);
    angles(ti,:) = trial(angleOpts(target,:));
end

channelCenters = 0:(2*pi/(N)):(2*pi-pi/(N));

%% Compute the likelihood of each trial
for tt = 1:length(fixedParams.trialTypes)
    
    idxs = adata(:,2)==fixedParams.trialTypeVals(tt);
    % get trials with this trial type
    tdata = adata(idxs,:);
    trial_likelihoods = nan(size(tdata,1),1);
    % for each trial, re-organize the angles so we can get the likelihoods
    % each column is the target, side, feat, and irrelevant angles all in
    % one calculation
    angleOpts = [8     9    10    11
                9     8    11    10
                10    11     8     9
                11    10     9     8];
    reOrder = angleOpts-7;

    angles = zeros(size(tdata,1),4);
    for ti = 1:size(tdata,1)
        angles(ti,:) = trial(angleOpts(tdata(ti,1),:));
    end

    % rotate all the angles relative to the response angle 
    angles = mod(angles-repmat(tdata(:,12),1,4),2*pi);
    % wtf was this code doing? Maybe you can do this with the original
    % model, but with this one it breaks the way the channels work... 
%     angles = angdist(repmat(tdata(:,12),1,4),angles);

    % get the parameters for this trialtype
    if fixedParams.one_sensitivity && fixedParams.shared_sensitivity
        dt = params.d_sh;
    elseif fixedParams.shared_sensitivity
        dt = params.dt_sh;
    else
        if fixedParams.cued_sensitivity
            if tt==1
                dt = params.dt_1;
            else
                dt = params.dt_cu;
            end
        else
            dt = params.(sprintf('dt_%i',tt));
        end
    end
    
    minmax = 5*dt;
    rL = 200;
    rrange = linspace(-minmax,1+minmax,rL);
    rrangec = rrange(:);
    
    dr = rrange(2)-rrange(1);

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
    
    %%
    % get the weights organized
    weights = [dt bs bf bi]';
    
    % pre-computed the repmat
    channelCentersRep = repmat(channelCenters',1,4);
    
    % run model
    clear cResp
    for ai = 1:size(angles,1)
        % for each trial, compute the channel responses
        % first we find the distances from all channels, and weight this by the
        % dprime vector
        cResp = (1-pscale_noconvert(angdist(channelCentersRep,repmat(angles(ai,:),100,1))*180/pi))*weights;
        % we now have the channel resp function, so the question is now how
        % likely was the real response to occur. The real response was rotated
        % to zero so we just have to compute the probability that the channel
        % at zero would have won over all others

        % we'll use rrange to do this, basically what we have to do is compute
        % the normpdf function across rrange for the encoder at zero, and then
        % multiply each of these values by the normcdf of all the other
        % channels for < that value
        % note that sigma = 1
        % get the probability range
        channel0pdf = normpdf(rrange,cResp(1),1)*dr;
        % compute the normcdf for all the other channels   
        allcdf = normcdf(repmat(rrangec,1,N-1),repmat(cResp(2:end)',rL,1),1)*dr;
        % now compute the probability by taking the prod across channels
        % for allcdf and multiplying by the channel0pdf values
        trial_likelihoods(ai) = channel0pdf * prod(allcdf,2);
        
        % test code: compute and plot the full likelihood
%         clear fulllike
%         for chan = 1:length(cResp)
%             channelpdf = normpdf(rrange,cResp(chan),1)*dr;
%             allcdf = normcdf(repmat(rrangec,1,N-1),repmat(cResp(setdiff(1:N,chan))',rL,1),1)*dr;
%             fulllike(chan) = channelpdf * prod(allcdf,2);
%         end
%         figure;
%         subplot(211);
%         plot(cResp);
%         subplot(212);
%         plot(fulllike);
    end
    %%
    probs(idxs) = trial_likelihoods;

    if any(trial_likelihoods==0)
        disp('Error');
        stop = 1;
    end
end

%%

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
    disp('todo');
    % compute the model output for this set of parameters, e.g. the 
%     fit.outs = outs;
%     fit.out = out;
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