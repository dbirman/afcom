function like = computeTCCfromPopulation(rads,sigma,p)
%% COMPUTE TCCPDF
% Usage: likelihood = computeTCCPDF(xs,dprime)
%
% Inputs:
%   rads - range of x values to evaluate at in *radians*
%   dprime - d prime value of the model
%
% Outputs:
%   likelihood - the likelihood that the TCC model will sample each x,
%                normalized

if iscolumn(rads)
    rads = rads';
end

npop = 100;
rpeaks = ((2*pi)/npop/2-pi):(2*pi/npop):pi;

responses = zeros(length(rads),length(rpeaks));
for ri = 1:length(rads)
    responses(ri,:) = normpdf(rpeaks,rads(ri),sigma);
end

% weight the responses by the exponential
weights = exp(-p*abs(rpeaks));

like = responses * weights';

%% test script
% x = -pi:pi/128:pi;
% 
% stim = 0;
% width = pi/6;
% 
% npop = 100;
% rpeaks = ((2*pi)/npop/2-pi):(2*pi/npop):pi;
% 
% response = normpdf(rpeaks,rads,width);
% 
% weights = exp(-2*abs(rpeaks));
% 
% % build the responses
% for ci = 0:(length(weights)-1)
%     out(ci+1) = response * circshift(weights,[ci,0])';
% end
% 
% figure(1);
% clf
% hold on
% % plot(rpeaks,response,'ok');
% out = response.*weights;
% plot(rpeaks,out./sum(out),'-r');
% 
% % now plot the TCC model prediction
% plot(rpeaks,computeTCCPDF(rpeaks,1.9),'-b');
% 
% legend({'Population readout','TCC'});