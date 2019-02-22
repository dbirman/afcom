function ac_plotSingleKappa(kappa)

h = figure;
    x = -pi:pi/32:pi;
    y = vonMises(x,0,kappa);
%% plot the vonmises
subplot(2,1,1);
plot(x,y,'-k');
axis([-pi pi 0 2]);
set(gca,'XTick',-pi:pi/2:pi,'XTickLabel',{'-180','-90','0','90','180'});
set(gca,'YTick',[0 1]');
axis square;
title(sprintf('FWHM (degs): %i',round(2*sqrt(2*log(2))*180/pi/sqrt(kappa))));
drawPublishAxis;

%% Plot a circle and add lines according to the distribution, put 0 north
r = 1;
off = pi/2;

subplot(2,1,2); hold on
circles(0,0,r,'facecolor','none');

plot([0 0],[0 1],'--r');
axis([-1 1 -1 1]*2.6);
axis square
for xi = 1:length(x)
    % get the position
    x1 = r * cos(x(xi)+off);
    y1 = r * sin(x(xi)+off);

    x2 = (r+y(xi))*cos(x(xi)+off);
    y2 = (r+y(xi))*sin(x(xi)+off);

    plot([x1 x2],[y1 y2],'-k','LineWidth',1);
end
drawPublishAxis('figSize=[10,7]');

savepdf(h,fullfile('~/proj/afcom/figures/kappas',sprintf('kappa_%if',round(kappa))));