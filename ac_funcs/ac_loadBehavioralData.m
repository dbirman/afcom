function [header,adata] = ac_loadBehavioralData(sid)
%%
files = dir(fullfile('~/data/afcom/',sid,'*.mat'));

maxTrackLength = 0;
for fi = 1:length(files)
    load(fullfile('~/data/afcom/',sid,files(fi).name));
    exp = getTaskParameters(myscreen,task);
    e{fi} = exp{1};
    mt{fi} = stimulus.data.mouseTrack(1:e{fi}.nTrials,:);
    mt{fi}(mt{fi}==0) = nan;
    maxTrackLength = max(maxTrackLength,size(mt{fi},2));
end

%% concatenate all trials
avars = {'target','trialType','cue','duration','dead','targetAngle','distractorAngle','angle1','angle2','angle3',...
    'angle4','respAngle','respDistance','distDistance'};

for ai = 1:length(avars)
    eval(sprintf('%s = [];',avars{ai}));
end

runs = [];
reactionTime = [];
runcount = [0 0];

for run = 1:length(e)
    runs = [runs run*ones(1,e{run}.nTrials)];
    runcount(e{run}.parameter.cue(1)) = runcount(e{run}.parameter.cue(1)) + 1;
    for ai = 1:length(avars)
        val = avars{ai};
        if isfield(e{run}.parameter,avars{ai})
            eval(sprintf('%s = [%s e{run}.parameter.%s];',val,val,val));
        else
            eval(sprintf('%s = [%s e{run}.randVars.%s];',val,val,val));
        end
    end
    reactionTime = [reactionTime e{run}.reactionTime];
end

%% create one giant matrix, but just of a few variables that matter
header = [avars {'runs' 'RT'}];
data = [];
for ai = 1:length(avars)
    data = eval(sprintf('[data %s'']',avars{ai}));
end
data = [data runs' reactionTime'];
adata = data;