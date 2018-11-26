function fit = ac_fitVonMises(adata,mode)
%% AC_FITVONMISES
% Fit a von mises to each condition with lapse rates
%
% MODE OPTIONS
%   'lapseall'; fits a lapse rate separately to each condition
global fixedParams

fixedParams.trialTypes = {'all','spatial','feature','target','baseline'};
fixedParams.trialTypeVals = [0 1 2 3 4];
fixedParams.lapseall = ~isempty(strfind(mode,'lapseall'));

%% PARAMETERS

for tt = 1:length(fixedParams.trialTypes)
    params.(sprintf('kappa%i',tt)) = [1 0 inf];
    if fixedParams.lapseall
        params.(sprintf('lapse%i',tt)) = [1 0 inf];
    end
end

if ~fixedParams.lapseall
    params.lapse = [1 0 inf];
end

% don't fit a bias
params.bias = 0; % [0 -pi pi];

[ip,minp,maxp] = initParams(params);

options = optimoptions('fmincon','Algorithm','active-set','TolFun',1,'TolCon',1,'Display','off'); % set a limit or it goes on foreeeeeeeeeeeever

bestparams = fmincon(@(p) vmlike(p,adata,0),ip,[],[],[],[],minp,maxp,[],options);
[~,fit] = vmlike(bestparams,adata,1);

function [likelihood, fit] = vmlike(params,adata,computeOutput)
global fixedParams

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
    
    if fixedParams.lapseall
        lapse = params.(sprintf('lapse%i',tt));
    else
        lapse = params.lapse;
    end
    
    kappa = params.(sprintf('kappa%i',tt));
    
    for ti = 1:size(tdata,1)
        trial = tdata(ti,:);

        if ~trial(5)
            % not dead
            probs(count) = lapse * lapseProb + (1-lapse) * vonMises(trial(12),trial(6),kappa);
            count = count + 1;
        end
    end
end

likelihood = -nansum(log(probs));

if isnan(likelihood)
    warning('Likelihood evaluated to NaN');
    stop = 1;
end

fit.trialProbs = probs;
fit.likelihood = likelihood;
fit.params = params;

if computeOutput
    
%     x = -pi:pi/128:pi;
%     fit.x = x;
%     fit.out = params.lapse * lapseProb + (1-params.lapse) * vonMises(x,0,params.kappa);
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