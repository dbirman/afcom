function [header,adata] = aca_loadBehavioralData()
%%
files = dir(fullfile('~/data/afcom_avg/',mglGetSID,'*.mat'));

maxTrackLength = 0;
for fi = 1:length(files)
    load(fullfile('~/data/afcom_avg/',mglGetSID,files(fi).name));
    exp = getTaskParameters(myscreen,task);
    e{fi} = exp{1};
    mt{fi} = stimulus.data.mouseTrack(1:e{fi}.nTrials,:);
    mt{fi}(mt{fi}==0) = nan;
    maxTrackLength = max(maxTrackLength,size(mt{fi},2));
end

% clear duration
% warning('adding duration = 1 if missing');
% for ei = 1:length(e)
%     if ~isfield(e{ei}.parameter,'duration')
%         e{ei}.parameter.duration = ones(size(e{ei}.parameter.trialType));
%     end
% end

%% concatenate all trials
pvars = {'target'};
rvars = {'dead','duration','trialType','targetAngle','respAngle','respDistance','target1','target2','angle1','angle2','angle3','angle4'};
runs = [];

header = [pvars rvars];

for pii = 1:length(pvars)
    eval(sprintf('%s = [];',pvars{pii}));
end
for ri = 1:length(rvars)
    eval(sprintf('%s = [];',rvars{ri}));
end

runcount = [0 0];
for run = 1:length(e)
    if e{run}.nTrials>0
        runs = [runs run*ones(1,e{run}.nTrials)];
        runcount(e{run}.parameter.cue(1)) = runcount(e{run}.parameter.cue(1)) + 1;
        for pii = 1:length(pvars)
            eval(sprintf('%s = [%s e{run}.parameter.%s];',pvars{pii},pvars{pii},pvars{pii}));
        end
        for ri = 1:length(rvars)
            eval(sprintf('%s = [%s e{run}.randVars.%s];',rvars{ri},rvars{ri},rvars{ri}));
        end
    end
end

eval('dur = duration;');

% copy the target info
for ai = 1:length(target1)
    eval(sprintf('target1(ai) = angle%i(ai);',target1(ai)));
    eval(sprintf('target2(ai) = angle%i(ai);',target2(ai)));
end

%% concatenate mouse tracks
% amt = nan(length(target),maxTrackLength);
% start = 1;
% for run = 1:length(e)
%     stop = (start+e{run}.nTrials-1);
%     amt(start:stop,1:size(mt{run},2)) = mt{run};
%     start = stop + 1;
% end
% 
% %% go backward through mouseTracks and fix jumps
% % assume that you end near zero, so if you jump -pi you need to -pi the
% % earlier section, etc
% amt = fliplr(amt);
% for ai = 1:size(amt,1)
%     track = amt(ai,:);
%     dtrack = diff(track);
%     posidx = find(dtrack>5);
%     negidx = find(dtrack<-5);
%     for pii = 1:length(posidx)
%         idx = posidx(pii)+1;
%         track(idx:end) = track(idx:end)-2*pi;
%     end
%     for nii = 1:length(negidx)
%         idx = negidx(nii)+1;
%         track(idx:end) = track(idx:end)+2*pi;
%     end
%     dtrack = diff(track);
%     amt(ai,:) = track;
% end
% amt = fliplr(amt);

%% create one giant matrix, but just of a few variables that matter

% pvars = {'target'};
% rvars = {'dead','duration','trialType','targetAngle','respAngle','respDistance'};

data = [runs' dead' trialType' respDistance' dur' target1' target2'];
keepIdxs = ~data(:,2);
data = data(keepIdxs,:);

adata = data;

disp(sprintf('Total trials: %i',size(data,1)));