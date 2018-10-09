%% AFCOM: Basic Von Mises Model
% This is a simple model for the afcom behavioral data. We assume that each
% stimulus is encoded as a Von Mises distribution and that the likelihood
% function is a weighted sum of these.
%
% For each of the five conditions in the experiment we fit 4 weights
% corresponding to the normalized weights on each stimulus:
%  (1) The target stimulus
%  (2) The same-side stimulus
%  (3) The feature-matched stimulus
%  (4) The off feature/side stimulus
% Along with the weights we fit the four K concentration parameters 
%
% Future to-do:
%
% Bias parameters can be added for (1) response bias to particular regions
% (encoded by cosine basis functions), (2) choice history bias

[headers,adata,amt] = ac_loadBehavioralData('s300');

%% Call the fitting routine

% first fit on the cue direction data
colordata = adata(adata(:,3)==1,:);
colorfit = ac_fitVM(colordata,headers);

dirdata = adata(adata(:,3)==2,:);
dirfit = ac_fitVM(dirdata,headers);

%% temp figure
fit = colorfit;

x = -pi:pi/128:pi;
    
    
h = figure;
plot(x,fit.out');
vline(fit.tAngle);
text(fit.tAngle,0.02,'Target');
vline(fit.sAngle);
text(fit.sAngle,0.02,'Same-side');
vline(fit.fAngle);
text(fit.fAngle,0.02,'Same-feature');
vline(fit.dAngle);
text(fit.dAngle,0.02,'Distractor');
legend({'All','Side','Feature','Target','Baseline'});
title('Color estimation');
drawPublishAxis;

savepdf(h,fullfile('~/proj/afcom/report_color.pdf'));

fit = dirfit;


h = figure;
plot(x,fit.out');
vline(fit.tAngle);
text(fit.tAngle,0.02,'Target');
vline(fit.sAngle);
text(fit.sAngle,0.02,'Same-side');
vline(fit.fAngle);
text(fit.fAngle,0.02,'Same-feature');
vline(fit.dAngle);
text(fit.dAngle,0.02,'Distractor');
legend({'All','Side','Feature','Target','Baseline'});
title('Direction estimation');
drawPublishAxis;
savepdf(h,fullfile('~/proj/afcom/report_dir.pdf'));
%% Evaluate fit
ac_plotVMFit(fit);

