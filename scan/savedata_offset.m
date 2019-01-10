function savedata_offset(cfolder,version)

% (1) Extract the left and right flat map voxel timecourses from the
% concatenation timeseries
% (2) Pull PRF analysis and R^2 values
% (3) Save everything into ~/Box Sync/ATTFIELD

warning('Not checking version numbers!');
%% Move to Data WD
mrQuit
cd(fullfile('~/data/offsetRetinotopy/',cfolder));
folders = dir(pwd);
skip = 1;
for fi = 1:length(folders)
    if ~isempty(strfind(folders(fi).name,'Concatenation')), skip = 0; end
end
if skip
    disp(sprintf('Data folder %s has not been prepared for analysis',cfolder));
    return
end

%% Setup a view + Load Concatenation
view = newView();
view = viewSet(view,'curGroup','Concatenation');

%% Load the pRF analysis
view = loadAnalysis(view,sprintf('pRFAnal/%s','pRF')); % check analysis name!
pRF = view.analyses{1};

%% Get ROI coordinates
ROIs = {'V1','V2','V3','V4','V3a','V3b','V7','MT'};
for ri = 1:length(ROIs)
    coords{ri} = getROICoordinates(view,ROIs{ri});
end

%% Load the tSeries
for scan = 1:3
    view = viewSet(view,'curScan',scan); % make sure scan # is correct
    
    tSeries{scan} = loadTSeries(view,view.curScan,[],[],[],[],[],view.curGroup);
end

%% Save
data.tSeries = tSeries;
data.coords = coords;
data.pRF = pRF;
data.ROIs = ROIs;
data.scan = 1:3;
data.group = view.curGroup;
data.cfolder = cfolder;
data.version = version;

save(fullfile('~/Box Sync/AFCOM_DATA',sprintf('data_%s.mat',cfolder)),'data','-v7.3');