function out = compTestChannelOut(units,stim,cK,mK)

% units is an array of unit properties

% stim is the stimulus properties

out = zeros(size(units,1),1); % this will be the "firing rate"

mult = stim(1)==units(:,1); % logical depending on whether a unit fires

coldist = angdist(units(:,2),stim(2));
motdist = angdist(units(:,3),stim(3));

% filter responses that don't fire at all
filt = (coldist>cK) .* (motdist>mK);

colresp = resp(coldist/cK);
motresp = resp(motdist/mK);

out = mult .* filt .* colresp .* motresp;

function r = resp(v)

if v>1
    r = 0;
else
    r = cos(v*pi/2).^4;
end