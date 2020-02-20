%% TCC ANALYSIS: AFCOM
% Load the data and fit the TCC model.

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
    calls = {'bads,spatial,feature,sh_sens,sh_bias',...
        'bads,spatial,feature,sh_sens,',...
        'bads,spatial,feature,sh_bias',...
        'bads,spatial,feature',...
        'bads,all,sh_sens,sh_bias',...
        'bads,baseline,sh_sens,sh_bias'};
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

%% Set permutation parameters
repeats = 100;

%% Do all the fits
disp(sprintf('Running fits for %i runs',length(infos)));
% sort by the number of trials in each one
len = zeros(size(infos));
for ii = 1:length(infos)
    len(ii) = size(infos{ii}.data,1);
end
[~,idxs] = sort(len,'descend');
infos = infos(idxs);
% run the actual fitting procedure
parfor ii = 1:length(infos)
    disp(infos{ii}.call);
    infos{ii}.fit = ac_fitTCCModel(infos{ii}.data,infos{ii}.call);
    infos{ii}.fit.data = infos{ii}.data;
    infos{ii}.fit.call = infos{ii}.call;
    infos{ii}.fit.dataType = infos{ii}.dataType;
    % now fit a number of permutation models as well, saving just the cv
    % values
    infos{ii}.fit.perm = zeros(1,repeats);
    info = infos{ii};
    for ri = 1:repeats
        info.data(:,12) = info.data(randperm(size(info.data,1)),12);
        % note this uses a faster fitting routine (has a max number of
        % function evals, so it won't run for ever)
        permFit = ac_fitTCCModel_perms(info.data,info.call);
        infos{ii}.fit.perm(ri) = permFit.cv.likelihood;
    end
end
disp('Fits complete');


%% Sort all the fits

allfits = {}; fits = {};
for si = [-1 1:length(subjects)]
    fit = {};
    for ci = 1:length(cues)
        for di = 1:length(mindurs)
            for ii = 1:length(infos)
                if infos{ii}.ci==ci && infos{ii}.di==di && infos{ii}.subj==si
                    % this is a fit that was done for this subject, for
                    % this condition, for this difficulty level
                    
                    % now we have to sort by the call, first pull out the
                    % baseline and cue4 fits
                    if ~isempty(strfind(infos{ii}.fit.call,'all'))
                        fit{ci,1} = infos{ii}.fit; % cue 4
                    elseif ~isempty(strfind(infos{ii}.fit.call,'baseline'))
                        fit{ci,2} = infos{ii}.fit; % baseline
                    elseif ~isempty(strfind(infos{ii}.fit.call,'sh_sens,sh_bias'))
                        fit{ci,3} = infos{ii}.fit; % shared model
                    elseif ~isempty(strfind(infos{ii}.fit.call,'sh_sens')) && isempty(strfind(infos{ii}.fit.call,'sh_bias'))
                        fit{ci,4} = infos{ii}.fit; % shared sensitivity
                    elseif ~isempty(strfind(infos{ii}.fit.call,'sh_bias')) && isempty(strfind(infos{ii}.fit.call,'sh_sens'))
                        fit{ci,5} = infos{ii}.fit; % shared bias
                    elseif ~isempty(strfind(infos{ii}.fit.call,'spatial,feature')) && isempty(strfind(infos{ii}.fit.call,'sh_bias')) && isempty(strfind(infos{ii}.fit.call,'sh_sens'))
                        fit{ci,6} = infos{ii}.fit; % no shared parameters
                    end
                end
            end
        end
    end
    if si>-1
        fits{si} = fit;
    end
end

%% Save the fits

% over-write with cross-validated fits
save(fullfile('~/proj/afcom/tcc_data.mat'),'fits');
