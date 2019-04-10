%% Plot a kappa distribution using lines
k = 0:.001:20;
f = 2*sqrt(2*log(2))*180/pi./sqrt(k);

fwhm = [Inf 180 135 90 45];
kappas = zeros(size(fwhm));
for fi = 1:length(fwhm)
    i = find(f<=fwhm(fi),1);
    kappas(fi) = k(i);
end


% kappas = [0 1 5 10 20];
h = figure;
for ki = 1:length(kappas)
    kappa = kappas(ki);
    x = -pi:pi/32:pi;
    y = vonMises(x,0,kappa);

    %% plot the vonmises
    subplot(2,length(kappas),ki);
    plot(x,y,'-k');
    axis([-pi pi 0 2]);
    set(gca,'XTick',-pi:pi/2:pi,'XTickLabel',{'-180','-90','0','90','180'});
    set(gca,'YTick',[0 1]');
    axis square;
    title(sprintf('FWHM: %i',round(2*sqrt(2*log(2))*180/pi/sqrt(kappa))));
    drawPublishAxis;
    
    %% Plot a circle and add lines according to the distribution, put 0 north
    r = 1;
    off = pi/2;

    subplot(2,length(kappas),length(kappas)+ki); hold on
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
    drawPublishAxis('figSize=[20,15]');
    
end

savepdf(h,fullfile('~/proj/afcom/figures/kappa_comparison.pdf'));