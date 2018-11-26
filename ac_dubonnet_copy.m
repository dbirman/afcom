%% AFCOM COPY FILES FROM DUBONNET

dfolder = '~/data/afcom/';
local = '~/data/afcom/';

subjects = [300 366 376 377];% 378];

call = sprintf('rsync -azP gru@dubonnet.stanford.edu:~/data/afcom ~/data');
system(call);

clear SIDs
for si = 1:length(subjects)
    SIDs{si} = sprintf('s%i',subjects(si));
end