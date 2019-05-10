
%%
subjects = [300 388 392 393];

clear SIDs
for si = 1:length(subjects)
    SIDs{si} = sprintf('s%i',subjects(si));
end