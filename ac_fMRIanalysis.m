%% fMRI analysis

path = fullfile('~/data/afcom/s037420190301/');
cd(path);

% What we need to do now is open the dataset, and pull every voxel. We'll
% then fit the channel tuning model. 

% First, load the different stimulus files and 

%% Setup a view + Load Concatenation
view = newView();
view = viewSet(view,'curGroup','Concatenation');
view = viewSet(view,'curScan',1); % make sure scan # is correct

%%
view = loadAnalysis(view,sprintf('erAnal/%s','all')); % check analysis name!

rois = loadROITSeries(view,allROIs,view.curScan,view.curGroup,'keepNAN=true');







%% HELPER CODE: FIX STIMULUS FILES 
% The stimulus files in the original task were missing some of the
% variables (trialType, target) 
count = 1;
cumTrials = 0;
fnum = 1;

files = dir(fullfile(path,'Etc','*.mat'));
for fi = 1:length(files)
    if isempty(strfind(files(fi).name,'original')) && isempty(strfind(files(fi).name,'fixed'))
%         disp('ALREADY DONE');
%         break
        
        fname = fullfile(path,'Etc',files(fi).name);
        load(fname);
        e = getTaskParameters(myscreen,task);
        e = e{1};
        
        b = stimulus.blocks{end};
        
        disp(sprintf('Number of trials in the block structure (+1 more than actually shown): %i',b.trial));
        disp(sprintf('Number of trials in the task structure (actually shown): %i',e.nTrials));
        
        % get the trialType and target fort hese
        tt = [];
        t = [];
        for ti = 1:e.nTrials
            tt = [tt b.trialType(count)];
            t = [t b.target(count)];
            count = count + 1;
        end
        
        % copy this info into the block
        for ti = 1:e.nTrials            
            task{1}{1}.block(ti).parameter.trialType = tt(ti);
            task{1}{1}.block(ti).parameter.target = t(ti);
        end
        
        task{1}{1}.parameter.trialType = [0 1 2];
        task{1}{1}.parameter.target = [1 2 3 4];
        
        disp(sprintf('I have used the data from %i trials so far',count));
%         
        save(fullfile(path,'Etc',sprintf('190301_stim%02.0f.mat',fnum)),'stimulus','task','myscreen');
        fnum = fnum+1;
    end
end