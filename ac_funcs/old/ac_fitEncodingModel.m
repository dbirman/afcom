function fit = ac_fitEncodingModel(adata,mode,fit)
%% AC_FITENCODINGMODEL
% The encoding model assumes that features are each encoded by independent
% von mises distributions. The likelihood of a response in any given
% condition is a weighted sum of these different distributions.
%
% PARAMETERS
%
%   lambda - lapse rate (lambda) and response rate (1-lambda)
%   beta_side - probability of sampling from the correct side
%   beta_feat - probability of sampling from the correct feature within the
%               correct side
%   beta_dist - probability of sampling from the correct feature but on the
%               wrong side
%
%   All of these von mises are multiplied according to their various
%   weights to compute the final likelihood distribution
%
% MODE OPTIONS
%
%   'bdist' - adds the beta_dist parameter, otherwise this is dropped
%   'bias' - adds a bias parameter
%
global fixedParams

fixedParams.trialTypes = {'all','spatial','feature','target','baseline'};
fixedParams.trialTypeVals = [0 1 2 3 4];
fixedParams.bdist = ~isempty(strfind(mode,'bdist'));

%% TEST FIT
if ~isempty(strfind(mode,'eval'))
    % evaluate fit of existing model
    [~,fit] = vmlike(fit.rawparams,adata,1);
    return
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
        
        trainfit = ac_fitEncodingModel(train,strcat(mode,',nocv'));
        testfit = ac_fitEncodingModel(test,strcat(mode,',eval'),trainfit);
        aprobs = [aprobs; testfit.probs];
        like(ri) = testfit.likelihood;
    end
    
    fit = ac_fitEncodingModel(adata,strcat(mode,',nocv'));
    fit.cv.probs = aprobs;
    fit.cv.likelihood = -nansum(log(aprobs));
    return
end

%% PARAMETERS

for tt = 1:length(fixedParams.trialTypes)
    if ~isempty(strfind(mode,'multikappa'))
        types = {'target','side','feat','dist'};
        for ti = 1:length(types)
            params.(sprintf('kappa%i_%s',tt,types{ti})) = [5 0.1 100 0.5 20];
        end
        fixedParams.multikappa = true;
    else
        params.(sprintf('kappa%i',tt)) = [5 0.1 100 0.5 20];
        fixedParams.multikappa = false;
    end
    params.(sprintf('lapse%i',tt)) = [0.1 0 1 0 0.7];
    params.(sprintf('beta_side%i',tt)) = [0.75 0 1 0.5 1];
    params.(sprintf('beta_feat%i',tt)) = [0.75 0 1 0.5 1];
    if fixedParams.bdist
        params.(sprintf('beta_dist%i',tt)) = [0.1 0 1 0 0.7];
    else
        params.beta_dist = 0.5;
    end
end

% don't fit a bias
if ~isempty(strfind(mode,'bias'))
    params.bias = [0 -pi pi -pi/8 pi/8];
else
    params.bias = 0; % [0 -pi pi];
end

[ip,minp,maxp,plb,pub] = initParams(params);

% BADS VERSION
if strfind(mode,'bads')
    bestparams = bads(@(p) vmlike(p,adata,0),ip,minp,maxp,plb,pub); %[],[],[],[],minp,maxp,[],[]);

% FMINCON VERSION
else
    options = optimoptions('fmincon','Algorithm','active-set','TolFun',1,'TolCon',1,'Display','off'); % set a limit or it goes on foreeeeeeeeeeeever

    bestparams = fmincon(@(p) vmlike(p,adata,0),ip,[],[],[],[],minp,maxp,[],options);
end

[~,fit] = vmlike(bestparams,adata,1);

fit.rawparams = bestparams;
fit.trialTypes = fixedParams.trialTypes;
fit.trialTypeVals = fixedParams.trialTypeVals;

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
count = 1;

% setup lapse probability distribution
lapseProb = vonMises(0,0,0);

% compute different trial types separately
for tt = 1:length(fixedParams.trialTypes)
    trialType = fixedParams.trialTypeVals(tt);
    tdata = sel(adata,2,trialType);
    
    if ~fixedParams.multikappa
        kappa = params.(sprintf('kappa%i',tt));
    end
    lapse = params.(sprintf('lapse%i',tt));
    beta_side = params.(sprintf('beta_side%i',tt));
    beta_feat = params.(sprintf('beta_feat%i',tt));
    if fixedParams.bdist
        beta_dist = params.(sprintf('beta_dist%i',tt));
    else
        beta_dist = params.beta_dist;
    end
    
    for ti = 1:size(tdata,1)
        trial = tdata(ti,:);

        if ~trial(5)
            % not dead
            
            % get the relevant angles
            switch trial(1)
                case 1
                    side = 2;
                    feat = 3;
                    dist = 4;
                case 2
                    side = 1;
                    feat = 4;
                    dist = 3;
                case 3
                    side = 4;
                    feat = 1;
                    dist = 2;
                case 4
                    side = 3;
                    feat = 2;
                    dist = 1;
            end
            
            aIdxs = [8 9 10 11];
            
            % compute the likelihood for the target side
            if ~fixedParams.multikappa
                likeTarget = vonMises(trial(12),trial(aIdxs(trial(1))),kappa);
                likeSide = vonMises(trial(12),trial(aIdxs(side)),kappa);
                likeFeat = vonMises(trial(12),trial(aIdxs(feat)),kappa);
                likeDist = vonMises(trial(12),trial(aIdxs(dist)),kappa);
            else
                likeTarget = vonMises(trial(12),trial(aIdxs(trial(1))),params.(sprintf('kappa%i_target',tt)));
                likeSide = vonMises(trial(12),trial(aIdxs(side)),params.(sprintf('kappa%i_side',tt)));
                likeFeat = vonMises(trial(12),trial(aIdxs(feat)),params.(sprintf('kappa%i_feat',tt)));
                likeDist = vonMises(trial(12),trial(aIdxs(dist)),params.(sprintf('kappa%i_dist',tt)));
            end
            
            likeSide = beta_feat * likeTarget + (1-beta_feat) * likeSide;
            likeOff = beta_dist * likeFeat + (1-beta_dist) * likeDist;
            
            likeTotal = beta_side * likeSide + (1-beta_side) * likeOff;
            
            probs(count) = lapse * lapseProb + (1-lapse) * likeTotal;
            count = count + 1;
        end
    end
end

likelihood = -nansum(log(probs));

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
    
    % set some arbitrary angles
    tTheta = 0;
    sTheta = -pi/2;
    fTheta = pi*1/3;
    dTheta = pi*2/3;
    
    for tt = 1:length(fixedParams.trialTypes)
        if ~fixedParams.multikappa
            kappa = params.(sprintf('kappa%i',tt));
        end
        lapse = params.(sprintf('lapse%i',tt));
        beta_side = params.(sprintf('beta_side%i',tt));
        beta_feat = params.(sprintf('beta_feat%i',tt));
        if fixedParams.bdist
            beta_dist = params.(sprintf('beta_dist%i',tt));
        else
            beta_dist = params.beta_dist;
        end
        
        % code from above
%         % compute the likelihood for the target side
        if ~fixedParams.multikappa
            likeTarget = vonMises(x,tTheta,kappa);
            likeSide = vonMises(x,sTheta,kappa);
            likeFeat = vonMises(x,fTheta,kappa);
            likeDist = vonMises(x,dTheta,kappa);
        else
                likeTarget = vonMises(x,tTheta,params.(sprintf('kappa%i_target',tt)));
                likeSide = vonMises(x,sTheta,params.(sprintf('kappa%i_side',tt)));
                likeFeat = vonMises(x,fTheta,params.(sprintf('kappa%i_feat',tt)));
                likeDist = vonMises(x,dTheta,params.(sprintf('kappa%i_dist',tt)));
        end
% 
        likeSide = beta_feat * likeTarget + (1-beta_feat) * likeSide;
        likeOff = (1-beta_dist) * likeFeat + beta_dist * likeDist;
% 
        likeTotal = beta_side * likeSide + (1-beta_side) * likeOff;
        
        fit.out(tt,:) = lapse * lapseProb + (1-lapse) * likeTotal;
    end
end















%% Helper routines

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