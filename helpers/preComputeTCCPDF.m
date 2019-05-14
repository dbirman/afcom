function probs = preComputeTCCPDF(rads,dprimes)
%% RETURN PRECOMPUTED *interpolated* TCCPDF
% This function is a very fast version of computeTCCPDF. The first time
% it's called it will load a global which has the TCCPDF likelihood
% functions pre-computed over a reasonably wide range. It then uses interp2
% to interpolate the likelihood across this space. 

if length(dprimes)==1
    repmat(dprimes,size(rads));
end

global tccLike

if ~isfield(tccLike,'now') || ((now-tccLike.now)>(1/24))
    % clear if it's been more than an hour
    tccLike = struct;
end

fname = fullfile('~/data/afcom_avg/','tcc.mat');
% either create from scratch or load if the file exists
if ~isfile(fname)
    disp('Rebuilding TCCPDF');
    tccLike = struct;
    tccLike.xs = 0:pi/256:pi;
    tccLike.dprimes = [0:.05:2 2:.1:6 6:.4:10 10:20 20:2:50 50:5:100];
    tccLike.dprimes = sort(unique(tccLike.dprimes));
    tccLike.like = zeros(length(tccLike.dprimes),length(tccLike.xs));
    disppercent(-1/length(tccLike.dprimes));
    for di = 1:length(tccLike.dprimes)
        tccLike.like(di,:) = computeTCCPDF(tccLike.xs,tccLike.dprimes(di));
        disppercent(di/length(tccLike.dprimes));
    end
    disppercent(inf);
    % plot to see that it worked
%         figure;
%         surf(tccLike.xs,tccLike.dprimes,tccLike.like);
    % save
    save(fname,'-struct','tccLike');
elseif isempty(fields(tccLike))
    tccLike = load(fname);
end

probs = interp2(tccLike.xs,tccLike.dprimes,tccLike.like,rads,dprimes);

if any(isnan(probs))
    warning('tccLike.like is not sufficiently large to encompass the data! Some data evaluated to NaN: breaking to keyboard');
    stop = 1;
    keyboard
end