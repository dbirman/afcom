
%%
subjects = [300 366 376 377 378];

clear SIDs
for si = 1:length(subjects)
    SIDs{si} = sprintf('s%i',subjects(si));
end