
%% Subjects
subjects = [89 91 300 366 376 377 378];
% removed 90, insufficient vespene gas (actually, didn't eye track well). 

clear SIDs
for si = 1:length(subjects)
    SIDs{si} = sprintf('s%03.0f',subjects(si));
end

allfemale = [1 1 1 0 0 0 1 1 0 1 1 0 1 0 1 0];
allyob = [1992 1998 1994 1985 1982 1990 1996 2000 1999 1996 1998 1990 1997 1996 1997 1998];
disp(sprintf('%i female %i male',sum(allfemale),sum(~allfemale)));
disp(sprintf('Average age %2.0f, range %2.0f - %2.0f',mean(2019-allyob),min(2019-allyob),max(2019-allyob)));

%% Setup colormap

cmap_ = colorblindmap/255;

clear cmap;
cmap(5,:) = [0.5 0.5 0.5];
cmap(4,:) = cmap_(7,:);
cmap(3,:) = cmap_(4,:);
cmap(2,:) = cmap_(8,:);
cmap(1,:) = cmap_(1,:);