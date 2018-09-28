%% AFCOM COPY FILES FROM DUBONNET

dfolder = '~/data/afcom/';
local = '~/data/afcom/';

subjects = [300];

for si = 1:length(subjects)
    call = sprintf('scp -r gru@dubonnet.stanford.edu:~/data/afcom ~/data');
    system(call);
end