%% TCC ANALYSIS: AFCOM
% Load the data and fit the TCC model.

warning('THIS SCRIPT ONLY RUNES THE ONE-SENSITIVITY MODEL');

% this was set up wrong originally. In fact we need to fit the following
% models:

% (1) Shared model (already fit), i.e. cued and uncued have separate
% parameters (you have to add together the likelihoods for these to be
% comparable to the following models:)
% (2) Uncued + cued share sensitivity parameters - todo
% (3) Uncued + cued share bias parameters - todo
% (4) Uncued + cued share all parameters  - todo 

% In the original version I fit a single d' parameter for all four dot
% patches, but that's stupid (although it works). In fact we want to see if
% the model "needs" the sensitivity or bias parameters to explain the effect of
% cueing. 

ac_setup;

cues = [1 2];

alldata = [];
adatas = {};
infos = {};

for si = 1:length(subjects)
    %% Load data
    [headers,adata] = ac_loadBehavioralData(SIDs{si});
    
% avars = {'target','trialType','cue','duration','dead','targetAngle','distractorAngle','angle1','angle2','angle3',...
%     'angle4','respAngle','respDistance','distDistance'};

% TARGET TRIALTYPE CUE DURATION DEAD TARGETANGLE DISTRACTORANGLE
%    1      2       3     4       5       6           7
% ANGLE1 ANGLE2 ANGLE3 ANGLE4 RESPANGLE RESPDISTANCE DISTDISTANCE RUNS RT
%    8      9     10     11      12          13            14      15  16
   
    %% Remove dead trials
    disp(sprintf('Dropping %i trials that included eye movements',sum(adata(:,5))));
    adata = adata(~adata(:,5),:);
    
    %% Remove any trials with duration > 0.3 (these were "easy" trials in the initial versions)
    disp(sprintf('Dropping %i trials from training or easy mode',sum(adata(:,4)>0.3)));
    adata = adata(adata(:,4)<=0.3,:);
    
    %% Drop the first 100 trials
%     adata = adata(151:end,:);
%     disp(sprintf('Dropping 150 trials to account for training time'));

    
    %% Save
    adatas{si} = adata;
    l(si) = size(adata,1);
    
    %% Concatenate w/ all data
    if size(adata,1)>250
        disp(sprintf('Adding data to alldata'));
        alldata = [alldata ; adata];
    end
    
    %% Display count
    disp(sprintf('Subject %s: %i trials',SIDs{si},size(adata,1)));
    
    %% Setup model fits
    colordata = sel(adata,3,1); % report color
    dirdata = sel(adata,3,2);
    
    reportType = {'color','direction'};
    calls = {'bads,all,spatial,feature,cued_sens,sh_bias',...
        'bads,all,spatial,feature,cued_bias,sh_sens',...
        'bads,all,spatial,feature,cued_sens,cued_bias',...
        'bads,all,spatial,feature,sh_sens,sh_bias'}; % all shared
    mindurs = [0.25];
    maxdurs = [0.3];
    durationType = {'hard'};
    for ci = 1:length(cues)
        for di = 1:length(mindurs)
            for li = 1:length(calls)
                data = sel(adata,3,cues(ci));
                data = fil(data,4,'>=',mindurs(di));
                data = fil(data,4,'<=',maxdurs(di));
                if size(data,1)>0
                    info = struct;
                    info.data = data;
                    info.call = calls{li};
                    info.ci = ci;
                    info.di = di;
                    info.subj = si;
                    info.dataType = sprintf('%s report %s',durationType{di},reportType{ci});
                    infos{end+1} = info;
        %             fit{ci,di} = ac_fitEncodingModel(data,'nocv,bdist');%'nocv');
                    disp(sprintf('%s: %i trials',info.dataType,size(data,1)));
                else
                    disp(sprintf('%s: skipping',sprintf('%s report %s',durationType{di},reportType{ci})));
                end
            end
        end
    end
end

%% Do an all subject fit
di = 1;
allinfos = {};
for ci = 1:2
    data = sel(alldata,3,cues(ci));
    data = fil(data,4,'>=',mindurs(di));
    data = fil(data,4,'<=',maxdurs(di));
    for li = 1:length(calls)
        info = struct;
        info.data = data;
        info.call = calls{li};
        info.ci = ci;
        info.subj = inf;
        info.dataType = sprintf('hard report %s',reportType{ci});
        allinfos{end+1} = info;
    end
end
for ii = 1:length(allinfos)
    disp(allinfos{ii}.call);
    allinfos{ii}.fit = ac_fitTCCModel(allinfos{ii}.data,allinfos{ii}.call);
    allinfos{ii}.fit.call = allinfos{ii}.call;
    allinfos{ii}.fit.dataType = allinfos{ii}.dataType;
end
save(fullfile('~/proj/afcom/tcc_all_sens.mat'),'allinfos');

%% nocv!

load(fullfile('~/proj/afcom/tcc_all_sens.mat'));
for ii = 1:length(allinfos)
    allinfos{ii}.call = strcat(allinfos{ii}.call,',nocv');
end
%% Now attempt to fit the "corrected" model

% parfor ii = 1:4%length(allinfos)
%     disp(allinfos{ii}.call);
%     allinfos{ii}.cfit = ac_fitTCCModel_corrected(allinfos{ii}.data,allinfos{ii}.call);
%     allinfos{ii}.cfit.call = allinfos{ii}.call;
%     allinfos{ii}.cfit.dataType = allinfos{ii}.dataType;
% end
% save(fullfile('~/proj/afcom/tcc_all_sens_corrected.mat'),'allinfos');

%% Do all the fits
disp(sprintf('Running fits for %i runs',length(infos)));
% sort by the number of trials in each one
len = zeros(size(infos));
for ii = 1:length(infos)
    len(ii) = size(infos{ii}.data,1);
end
[~,idxs] = sort(len,'descend');
infos = infos(idxs);

%% run the actual fitting procedure
parfor ii = 1:length(infos)
    disp(infos{ii}.call);
    infos{ii}.fit = ac_fitTCCModel(infos{ii}.data,infos{ii}.call);
    infos{ii}.fit.call = infos{ii}.call;
    infos{ii}.fit.dataType = infos{ii}.dataType;
end
disp('Fits complete');
save(fullfile('~/proj/afcom/tcc_one_sens.mat'),'infos');

%%
% Save the new fits into the existing fit dataset
load(fullfile('~/proj/afcom/tcc_one_sens.mat'));
load(fullfile('~/proj/afcom/tcc_data.mat'));

% Re-organize infos into cued_fits

cued_fits = {};

for ii = 1:length(infos)
    info = infos{ii};
    
    if ~isempty(strfind(info.call,'cued_sens,cued_bias'))
        model = 4; % all cued separate from uncued
    elseif ~isempty(strfind(info.call,'sh_sens,sh_bias'))
        model = 1; % all shared
    elseif ~isempty(strfind(info.call,'cued_sens,sh_bias'))
        model = 2;
    else
        model = 3;
    end
    
    cued_fits{info.subj}{info.ci,model} = info.fit;
end

save(fullfile('~/proj/afcom/tcc_data.mat'),'fits','cued_fits');

%% Model recovery code


