function fit = aca_fitTCCPopulationModel(adata,mode,fit)
%% ACA_FITTCCMODEL
% Note this is aca not ac, so it's fitting the perceptual averaging data
% (afcom_avg = aca)
%
% This version of the model fits the "TCC" model from this paper:
% Schurgin, M. W., Wixted, J. T., & Brady, T. F. (2018). Psychological scaling reveals a single parameter framework for visual working memory. bioRxiv, 325472.
%
%
% PARAMETERS
% 
%
%
%
% MODE OPTIONS
%

global fixedParams

fixedParams.trialTypes = {'spatial','feature'};
fixedParams.trialTypeVals = [1 2];

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

params.sigma = [pi/8 0 10*pi pi/16 pi];
params.p = [2 0 10 1 3];
% params.int = [1 0 10 0.25 1.75]; 
% params.typeInt = [0 -1 1 -0.4 0.4]; % this is the difference between feature and spatial (feature==1, spatial==0)*typeInt

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
fit.adata = adata;

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

% get the predicted likelihood for each response relative to the target
% angle
probs = preComputeTCCPDF(adata(:,4),params.p);
% probs = computeTCCfromPopulation(adata(:,4),params.sigma,params.p);


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