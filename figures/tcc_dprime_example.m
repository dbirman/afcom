h = figure;
hold on

x = 0:pi/32:pi;

dprimes = logspace(-1,0.5,6);
dprimes = fliplr(dprimes);
dprimes = dprimes([1 3 5]);

cmap = brewermap(13,'Reds');

dp = {};
for i = 1:length(dprimes)
    disp(dprimes(i));
    plot(pscale(x*180/pi),computeTCCPDF(x,dprimes(i)),'-','Color',cmap(14-i,:));
    dp{end+1} = sprintf('d'' = %1.1f',round(dprimes(i)*10)/10);
end
xlabel('Response distance from target (normalized psychophysical distance');
ylabel('Probability of response (pdf)');
axis([0 1 0 0.3]);
set(gca,'XTick',0:.2:1,'YTick',0:.1:.3);
legend(dp);

drawPublishAxis('figSize=[25,10]','poster=1');

savepdf(h,fullfile('~/proj/afcom/figures','dprimes_examples.pdf'));