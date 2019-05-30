%% Simulate data from the TCC model and explore whether the model routines can recover the parameters

nTrials = 3000;

dpt = 2; % dprime of the target dot patch
dps = 0.5; % dprime of the same-side
dpf = 0.5; % dprime of the same-feature
dpi = 0.5; % dprime of the irrelevant

dprimes = [dpt dps dpf dpi];

bs = 0.9; % probability of being on the correct side
bf = 0.9; % probability of getting the feature right | on the correct side
bi = 0.5; % probability of getting the feature wrong | on the wrong side

probs = [bs*bf bs*(1-bf) (1-bs)*(1-bi) (1-bs)*bi];

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

% now compute the full likelihood by collapsing over the four probabilities
% for each trial, again we'll do it is a loop for clarity
liket = zeros(size(likes,1),length(xs));
for ai = 1:size(likes,1)
    liket(ai,:) = squeeze(likes(ai,:,:))' * probs';
end
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
for ai = 1:size(liket,1)
    cs = cumsum(liket(ai,:));
    respAngle = interp1(cs,xs,rand*max(cs));
    adata(ai,:) = [1 0 0 1 0 angs(ai,1) nan angs(ai,:) respAngle angdist(respAngle,angs(ai,1)) nan 1 1];
end

%% Fit the data with the model:
fit = ac_fitTCCModel(adata,'nocv,bads,recovery');

%% Check how well we recovered the parameters
fit.params