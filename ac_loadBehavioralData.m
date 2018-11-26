function [header,adata, amt] = ac_loadBehavioralData(sid)
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
pvars = {'target','trialType','cue','duration'};
rvars = {'dead','targetAngle','distractorAngle','angle1','angle2','angle3',...
    'angle4','respAngle','respDistance','distDistance'};

for pii = 1:length(pvars)
    eval(sprintf('%s = [];',pvars{pii}));
end
for ri = 1:length(rvars)
    eval(sprintf('%s = [];',rvars{ri}));
end
runs = [];
reactionTime = [];

runcount = [0 0];
for run = 1:length(e)
    runs = [runs run*ones(1,e{run}.nTrials)];
    runcount(e{run}.parameter.cue(1)) = runcount(e{run}.parameter.cue(1)) + 1;
    for pii = 1:length(pvars)
        eval(sprintf('%s = [%s e{run}.parameter.%s];',pvars{pii},pvars{pii},pvars{pii}));
    end
    for ri = 1:length(rvars)
        eval(sprintf('%s = [%s e{run}.randVars.%s];',rvars{ri},rvars{ri},rvars{ri}));
    end
    reactionTime = [reactionTime e{run}.reactionTime];
end

%% concatenate mouse tracks
amt = nan(length(target),maxTrackLength);
start = 1;
for run = 1:length(e)
    stop = (start+e{run}.nTrials-1);
    amt(start:stop,1:size(mt{run},2)) = mt{run};
    start = stop + 1;
end

%% go backward through mouseTracks and fix jumps
% assume that you end near zero, so if you jump -pi you need to -pi the
% earlier section, etc
amt = fliplr(amt);
for ai = 1:size(amt,1)
    track = amt(ai,:);
    dtrack = diff(track);
    posidx = find(dtrack>5);
    negidx = find(dtrack<-5);
    for pii = 1:length(posidx)
        idx = posidx(pii)+1;
        track(idx:end) = track(idx:end)-2*pi;
    end
    for nii = 1:length(negidx)
        idx = negidx(nii)+1;
        track(idx:end) = track(idx:end)+2*pi;
    end
    dtrack = diff(track);
    amt(ai,:) = track;
end
amt = fliplr(amt);

%% create one giant matrix, but just of a few variables that matter
header = [pvars rvars {'runs' 'RT'}];
data = [];
for pii = 1:length(pvars)
    data = eval(sprintf('[data %s'']',pvars{pii}));
end
for ri = 1:length(rvars)
    data = eval(sprintf('[data %s'']',rvars{ri}));
end
data = [data runs' reactionTime'];
adata = data;
% keepIdxs = ~any(isnan(data(:,4)),2);
% adata = data(keepIdxs,:);
% amt = amt(keepIdxs,:);