ls%% AFCOM COPY FILES FROM DUBONNET

dfolder = '~/data/afcom/';
local = '~/data/afcom/';

call = sprintf('rsync -azP gru@dubonnet.stanford.edu:~/data/afcom ~/data');
system(call);