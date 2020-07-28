%% Load
load(fullfile('~/proj/afcom/tcc_all_sens.mat'));

clear allfits
for ii = 1:length(allinfos)
    info = allinfos{ii};
    
    if ~isempty(strfind(info.call,'cued_sens,cued_bias'))
        model = 4; % all cued separate from uncued
    elseif ~isempty(strfind(info.call,'sh_sens,sh_bias'))
        model = 1; % all shared
    elseif ~isempty(strfind(info.call,'cued_sens,sh_bias'))
        model = 2;
    else
        model = 3;
    end
    
    allfits{info.ci,model} = info.fit;
end

%%
% get the likelihoods
clear likes
for ci = 1:2
    for mi = 1:4
        likes(ci,mi) = allfits{ci,mi}.cv.likelihood;
    end
end

% MODEL 1 = ALL SHARED
% MODEL 2 = SEPARATE CUED SENSITIVITY
% MODEL 3 = SEPARATE CUED BIAS
% MODEL 4 = BOTH CUED SEPARATE
likes(1,2:4) - likes(1,1)
likes(2,2:4) - likes(2,1)