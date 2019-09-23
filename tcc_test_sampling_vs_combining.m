% This script is for checking that "sampling" from a set of likelihood
% distributions is equivalent to combining them in a weighted manner. It's
% possible these aren't equivalent since the distributions can have
% different peak values? 

%% Generate a set of likelihood functions using preComputeTCCPDF

nTrials = 5;

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

% generate the four likelihood functions

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

% size of likes is (ntrials | dot patch | likelihood)

% now either combine the distributinos by multiplying:
like_mult = zeros(size(likes,1),length(xs));
for ai = 1:size(likes,1)
    like_mult(ai,:) = squeeze(likes(ai,:,:))' * probs';
end

% and then compute choicesadata = zeros(0,16);
acs = cumsum(liket,2);
for ai = 1:size(liket,1)
    cs = acs(ai,:);
    for i = 1:100
        respAngle(ai) = interp1(cs,xs,rand*max(cs));
        
        temp(i) = respAngle;
    end
end

% or compute it 