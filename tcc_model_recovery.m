function [orig_params, fit_params, diff_params] = tcc_model_recovery()
%% Simulate data from the TCC model and explore whether the model routines can recover the parameters

nTrials = 500;

% dpt = 2; % dprime of the target dot patch
% dps = 1; % dprime of the same-side
% dpf = 1; % dprime of the same-feature
% dpi = rand*2+1; % dprime of the irrelevant
dpt = 2;
dps = 0.5;
dpf = 0.5;
dpi = 0.5;

dprimes = [dpt dps dpf dpi];
disp(dprimes);

% bs = rand*0.5+0.5; % probability of being on the correct side
% bf = rand*0.5+0.5; % probability of getting the feature right | on the correct side
% bi = rand*0.2; % probability of getting the feature wrong | on the wrong side
bs = 0.95;
bf = 0.9;
bi = 0.5;

probs = [bs*bf bs*(1-bf) (1-bs)*(1-bi) (1-bs)*bi];
disp(probs);

orig_params = [dprimes bs bf bi];
% therefore:
% bs*bf = probability of choosing target
% bs*(1-bf) = probability of choosing same-side
% (1-bs)*(1-bi) = probability of choosing same-feature
% (1-bs)*bi = probability of choosing the irrelevant distractor

%% Simulate trials

% generate four random angles for each trial
angs = rand(nTrials,4)*pi;

% generate the range of possible errors
xs = 0:pi/64:pi;

likes = zeros(size(angs,1),size(angs,2),length(xs));

% what we need to do now is get the full likelihood distribution (over the
% range 0:pi) for each of the angles presented. We'll sample from this
% distribution to create trials. For clarity, we'll do this as a for loop.

% compute the likelihood of each angle, given the dprime values. Note that
% because the TCC likelihood functions are slow to build this function
% pre-computes them and then interpolates.
for di = 1:size(angs,2)
    for ai = 1:size(angs,1)
        likes(ai,di,:) = preComputeTCCPDF(angdist(xs-angs(ai,di),0),dprimes(di));
    end
end

figure;
subplot(221);
imagesc(squeeze(likes(:,1,:)));
title('Target');
subplot(222);
imagesc(squeeze(likes(:,2,:)));
title('Side');
subplot(223);
imagesc(squeeze(likes(:,3,:)));
title('Feature');
subplot(224);
imagesc(squeeze(likes(:,4,:)));
title('Distractor');

% now compute the full likelihood by collapsing over the four probabilities
% for each trial, again we'll do it is a loop for clarity
liket = zeros(size(likes,1),length(xs));
for ai = 1:size(likes,1)
    liket(ai,:) = squeeze(likes(ai,:,:))' * probs';
end

figure;
imagesc(liket);
title('Full likelihood for each trial');

% sample uniformly from each trial to generate data
%   Columns 1 through 6
% 
%     'target'    'trialType'    'cue'    'duration'    'dead'    'targetAngle'
% 
%   Columns 7 through 11
% 
%     'distractorAngle'    'angle1'    'angle2'    'angle3'    'angle4'
% 
%   Columns 12 through 16
% 
%     'respAngle'    'respDistance'    'distDistance'    'runs'    'RT'
% 
% note that many of the columns aren't actually used by the fitting
% routines. The important stuff is target, trialType, cue, duration, dead,
% angle1-4, respAngle. 
adata = zeros(0,16);
acs = cumsum(liket,2);
for ai = 1:size(liket,1)
    cs = acs(ai,:);
    for i = 1:100
        respAngle = interp1(cs,xs,rand*max(cs));
        
        temp(i) = respAngle;
    end
    adata(ai,:) = [1 0 0 1 0 angs(ai,1) nan angs(ai,:) respAngle angdist(respAngle,angs(ai,1)) nan 1 1];
end

figure;
hist(adata(:,12));

%% Fit the data with the model:
fit = ac_fitTCCModel(adata,'nocv,bads,recovery');

%% Check how well we recovered the parameters

fit_params = [fit.params.dt_1 fit.params.ds_1 fit.params.df_1 fit.params.di_1 fit.params.bs_1 fit.params.bs_1 fit.params.bi_1];

diff_params = fit_params-orig_params;


return



%% Test script
N = 50;
orig = zeros(N,7);
fit = orig;
diff = orig;
parfor i = 1:N
    [orig(i,:) fit(i,:) diff(i,:)] = tcc_model_recovery;
    disp(i);
end
