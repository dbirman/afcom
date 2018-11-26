function ac_plotVMFit(fit,type,fname)
x = -pi:pi/128:pi;

cmap = brewermap(size(fit.out,1),'Dark2');
h = figure; hold on
for i = 1:size(fit.out,1)
    plot(x,fit.out(i,:)','Color',cmap(i,:));
end
vline(fit.tAngle);
text(fit.tAngle,0.02,'Target');
vline(fit.sAngle);
text(fit.sAngle,0.02,'Same-side');
vline(fit.fAngle);
text(fit.fAngle,0.02,'Same-feature');
vline(fit.dAngle);
text(fit.dAngle,0.02,'Distractor');
legend({'All','Side','Feature','Target','Baseline'});
title(sprintf('%s estimation',type));
axis([min(x) max(x) 0 0.025]);
drawPublishAxis;
savepdf(h,fullfile('~/proj/afcom/',sprintf('%s.pdf',fname)));