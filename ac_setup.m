
%% Subjects
subjects = [89 90 300 366 376 377 378];

clear SIDs
for si = 1:length(subjects)
    SIDs{si} = sprintf('s%03.0f',subjects(si));
end


%% Setup colormap

cmap_ = colorblindmap/255;

clear cmap;
cmap(5,:) = [0.5 0.5 0.5];
cmap(4,:) = cmap_(7,:);
cmap(3,:) = cmap_(4,:);
cmap(2,:) = cmap_(8,:);
cmap(1,:) = cmap_(1,:);