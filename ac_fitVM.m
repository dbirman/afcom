function fit = ac_fitVM(adata,header)
%% AFCOM: Fit VonMises model

global fixedParams

%% Data check

%% Setup parameters

% each task condition gets a kappa parameter and a beta weight

% target trials only have a single kappa value
params.kappa0 = [2 0 inf];
params.beta0 = 1;

% all the ohter conditions (1-4) have a beta weight and kappa for every
% target. We'll call them target, featdist, sidedist, distractor to keep
% track.
types = {'target','featdist','sidedist','distdist'};
for tt = 1:4
    for ti = 1:4
        type = types{ti};
        params.(sprintf('kappa_%i_%s',tt,type)) = [2 0 inf];
        if ti==1 % if target
            binit = 1; % set the target weight to 1 so that normalizing is simpler
        else
            binit = [0 -inf inf];
        end
        params.(sprintf('beta_%i_%s',tt,type)) = binit;
    end
end

%% Fit
[ip,minp,maxp] = initParams(params);

options = optimoptions('fmincon','Algorithm','active-set','TolFun',1,'TolCon',1,'Display','off'); % set a limit or it goes on foreeeeeeeeeeeever

bestparams = fmincon(@(p) vmlike(p,adata),ip,[],[],[],[],minp,maxp,[],options);
[finalLike,fit] = vmlike(bestparams,adata);

fit.params = getParams(bestparams);
fit.like = finalLike;

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
    
    switch trial(2)
        case 0
            % BASELINE TRIAL
            probs(ai) = vonMises(trial(11),trial(5),params.kappa0);
        otherwise
            % TARGET TRIAL
            % compute individial probabilities and weight
            probs(ai) = getWeightedProb(trial,params);
    end
end

likelihood = -sum(log(probs));

fit = struct;
% x = -pi:pi/128:pi;
% for i = 0:4
%     fit.out(i+1,:) = vonMises(x,0,params.(sprintf('kappa%i',i)));
% end


function prob = getWeightedProb(trial,params)

sidedists = [2 1 4 3];
featdists = [3 4 1 2];
distdists = [4 3 2 1];

sidedist = sidedists(trial(1));
featdist = featdists(trial(1));
distdist = distdists(trial(1));

% get the angle for each one
idx = [7 8 9 10];

targetAngle = trial(idx(trial(1)));
featAngle = trial(idx(featdist));
sideAngle = trial(idx(sidedist));
distAngle = trial(idx(distdist));

% get the probabilities
targetProb = vonMises(targetAngle,0,params.(sprintf('kappa_%i_target',trial(2))));
featProb = vonMises(targetAngle,0,params.(sprintf('kappa_%i_featdist',trial(2))));
sideProb = vonMises(targetAngle,0,params.(sprintf('kappa_%i_sidedist',trial(2))));
distProb = vonMises(targetAngle,0,params.(sprintf('kappa_%i_distdist',trial(2))));
probs = [targetProb featProb sideProb distProb];

% get the betas
betas = [params.(sprintf('beta_%i_target',trial(2))) ...
    params.(sprintf('beta_%i_featdist',trial(2))) ...
    params.(sprintf('beta_%i_sidedist',trial(2))) ...
    params.(sprintf('beta_%i_distdist',trial(2)))];
betas = betas ./ sum(betas);

prob = probs * betas';

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