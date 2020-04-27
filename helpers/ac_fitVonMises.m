% fit the TCC model to a single condition ONLY FOR THE AVERAGING TASK

% this is just a test function
function kappa = ac_fitVonMises(rads)

if any(isnan(rads))
    disp(sprintf('Dropping %i/%i NaNs',sum(isnan(rads)),length(rads)));
    rads = rads(~isnan(rads));
end
    
kappa = bads(@(x) -sum(log(evalVonMises(rads,x))),1,0,100,0.25,5);

function likes = evalVonMises(rads,k)

xs = 0:pi/128:pi;

% compute the PDF of a vonmises with params p
vm_pdf = vonMises(x,0,k);

likes = interp1(xs,vm_pdf,rads);

