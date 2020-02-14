xs = pi/64:pi/32:pi;

reportType = {'color','direction'};
resp = zeros(length(fits),2,5,length(xs));
model = resp;
for subj = 1:length(fits)
    for cond = 1:2
        sfits = fits{subj};
        cfits = {sfits{cond,1}, sfits{cond,2}, sfits{cond,3}};
        tts = {[0],[4],[1 2]};
        for ci = 1:length(tts)
            cfit = cfits{ci};
            tt = tts{ci};
            
            % get the adata from these fits, they're all the same so it
            % doesn't matter which one we use
            ttdata = [];
            for ti = 1:length(tt)
                ttdata = [ttdata ; sel(cfit.data,2,tt(ti))];
            end
        end
        
        % bin the values
        [c,~] = hist(ttdata(:,13),xs);
        % normalize counts (so we can average them later)
        c = c ./ sum(c);
        % save
        resp(subj,cond,tt+1,:) = c;
                 
            
        % compute the model fit
        dt = cfit.params.dt_sh;
        probs = computeTCCPDF(xs,dt);
        model(subj,cond,tt+1,:) = probs;
    end
end