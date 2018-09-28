function fit = ac_fitVM(adata,header)
%% AFCOM: Fit VonMises model

global fixedParams

%% Data check

%% Setup parameters

% each task condition gets a kappa parameter and a beta weight

% start by fitting only the target trials
params.kappa0 = [2 0 inf];
params.kappa1 = [2 0 inf];
params.kappa2 = [2 0 inf];
params.kappa3 = [2 0 inf];
params.kappa4 = [2 0 inf];
params.beta = 1;

%% Fit
[ip,minp,maxp] = initParams(params);

options = optimoptions('fmincon','Algorithm','active-set','TolFun',1,'TolCon',1,'Display','off'); % set a limit or it goes on foreeeeeeeeeeeever

bestparams = fmincon(@(p) vmlike(p,adata),ip,[],[],[],[],minp,maxp,[],options);
[finalLike,fit] = vmlike(bestparams,adata);

stop = 2;

function [likelihood,fit] = vmlike(params,adata)
%% Likelihood function
    
%   Columns 1 through 5
% 
%     'target'    'trialType'    'cue'    'dead'    'targetAngle'
% 
%   Columns 6 through 10
% 
%     'distractorAngle'    'angle1'    'angle2'    'angle3'    'angle4'
% 
%   Columns 11 through 13
% 
%     'respAngle'    'respDistance'    'distDistance'

params = getParams(params);


% get the predicted likelihood for each response relative to the target
% angle
for ai = 1:size(adata,1)
    trial = adata(ai,:);
    
    probs(ai) = vonMises(trial(11),trial(5),params.(sprintf('kappa%i',trial(2))));
end

likelihood = -sum(log(probs));

x = -pi:pi/128:pi;
for i = 0:4
    fit.out(i+1,:) = vonMises(x,0,params.(sprintf('kappa%i',i)));
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