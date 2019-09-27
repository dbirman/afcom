% fit the TCC model to a single condition

% this is just a test function
function dprime = fitTCC(rads)

if any(isnan(rads))
    disp(sprintf('Dropping %i/%i NaNs',sum(isnan(rads)),length(rads)));
    rads = rads(~isnan(rads));
end
    
dprime = bads(@(x) -sum(log(evalTCC(rads,x))),1,0,100,0.25,5);

function likes = evalTCC(rads,p)

xs = 0:pi/128:pi;
like = computeTCCPDF(xs,p);
likes = interp1(xs,like,rads);