%% Plot a TCC distribution using lines
dprimes = [0.6 0.95];%.25:2;


h = figure;
for di = 1:length(dprimes)
    dprime = dprimes(di);
    x = -pi:pi/32:pi;
    y = computeTCCPDF(x,dprime);

    %% plot the likelihood
    subplot(2,length(dprimes),di);
    plot(x,y,'-k');
    axis([-pi pi 0 0.5]);
    set(gca,'XTick',-pi:pi/2:pi,'XTickLabel',{'-180','-90','0','90','180'});
    set(gca,'YTick',[0 0.1]');
    axis square;
    title(sprintf('d'': %1.2f',dprime));
    drawPublishAxis;
    
    %% Plot a circle and add lines according to the distribution, put 0 north
    r = 3;
    off = pi/2;
    subplot(2,length(dprimes),length(dprimes)+di); hold on

    circles(0,0,1,'facecolor','none');
    
    y = y * 7.5;

    plot([0 0],[0 1],'--r');
    axis([-1 1 -1 1]*3);
    axis square
    for xi = 1:length(x)
        % get the position
        x1 = 1 * cos(x(xi)+off);
        y1 = 1 * sin(x(xi)+off);

        x2 = (1+y(xi))*cos(x(xi)+off);
        y2 = (1+y(xi))*sin(x(xi)+off);

        plot([x1 x2],[y1 y2],'-k','LineWidth',1);
    end
    drawPublishAxis('figSize=[20,15]');
    
end

savepdf(h,fullfile('~/proj/afcom/figures/dprime_comparison_duration.pdf'));