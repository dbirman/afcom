%% Simulate the TCC model, then fit the simulations
% The TCC model is from Schurgin et al. it's hard to implement because it
% requires a (deterministic) step that involves simulating the responses of
% abstract "units" in the brain. My implementation is even more problematic
% because I then go and do probabilistic weighting according to various
% possible sources of bias in subject responses. To validate the model this
% simulation code generates artificial datasets of a similar size to real
% datasets from observers and then asks whether the original parameters are
% recoverable.
%
% The main issue is that the probability parameters (bt/bs/bf) may trade
% off with the sensitivity parameters (dt/ds/df/di), which would be
% disasterous since it would mean the model is mis-specified. 