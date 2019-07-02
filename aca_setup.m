
%%
subjects = [82 87 300 388 392 393 409];

clear SIDs
for si = 1:length(subjects)
    SIDs{si} = sprintf('s03.0f',subjects(si));
end