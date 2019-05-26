%% Simulate data from the TCC model and explore whether the model routines can recover the parameters

nTrials = 1000;

dpt = 1; % dprime of the target dot patch
dps = 1; % dprime of the same-side
dpf = 1; % dprime of the same-feature
dpi = 1; % dprime of the irrelevant

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

xs = 0:pi/64:pi;
likes = zeros(size(angs,1),size(angs,2),length(xs));

% what we need to do now is get the full likelihood distribution (over the
% range 0:pi) for each of the angles presented. We'll sample from this
% distribution to create trials. 

% for clarity, we'll do this as a for loop.

% compute the likelihood of each angle, given the dprime values
for di = 1:size(angs,2)
    for ai = 1:size(angs,1)
        likes(ai,di,:) = preComputeTCCPDF(angdist(xs-angs(ai,di),0),dprimes(di));
    end
end

% now compute the full likelihood by collapsing over the four probabilities
% for each trial
liket = zeros(size(likes,1),length(xs));
for ai = 1:size(likes,1)
    liket(ai,:) = squeeze(likes(ai,:,:)) * probs';

% take the maximum 