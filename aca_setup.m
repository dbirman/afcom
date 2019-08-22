
%%
subjects = [82 87 300 388 392 393 409];

female = [1 1 0 1 0 1 0];
yob = [1996 1998 1990 1997 1996 1997 1998];

clear SIDs
for si = 1:length(subjects)
    SIDs{si} = sprintf('s03.0f',subjects(si));
end