ls%% AFCOM COPY FILES FROM DUBONNET

call = sprintf('rsync -azP gru@dubonnet.stanford.edu:~/data/afcom ~/data');
system(call);

dfolder = '~/data/afcom_avg/';

call = sprintf('rsync -azP gru@dubonnet.stanford.edu:%s ~/data',dfolder);
system(call);