function aca_plotTCCModel(fit)

adata = fit.adata;

%% plot the duration plot
h = figure; hold on
cmap = colorblindmap/255;

durs = 0.25:.05:.75;
dp = fit.params.int + fit.params.bdur*durs + fit.params.bdist*0.75*pi*0.5;

plot(durs,dp,'Color',cmap(2,:));
plot(durs,dp+fit.params.typeInt,'Color',cmap(3,:));

legend({'Cue side','Cue color'});
set(gca,'XTick',[0 0.25 0.5 0.75 1]);
xlabel('Duration (s)');
ylabel('dprime (a.u.)');
axis([0 1 0 1.5]);
title('Effect of duration at dist=67.5 degs');

drawPublishAxis('figSize=[10,5]');
% xs = 0:pi/64:pi;
% 
% dbins = linspace(0.25,0.75,8);
% for di = 2:length(dbins)    
%     low = dbins(di-1);
%     high = dbins(di);
%     mid = mean([low high]);
%     idxs = (adata(:,5)>=low).*(adata(:,5)<high);
%     dat = adata(logical(idxs),:);
%     
%     rads = dat(:,4);
%     
% end

savepdf(h,fullfile('~/proj/afcom/figures/duration.pdf'));

%% Plot the effect of distance
h = figure; hold on
cmap = colorblindmap/255;

dist = 0:pi/32:(0.75*pi);
dp = fit.params.int + fit.params.bdur*0.5 + fit.params.bdist*dist;

plot(dist,dp,'Color',cmap(2,:));
plot(dist,dp+fit.params.typeInt,'Color',cmap(3,:));

legend({'Cue side','Cue color'});
xlabel('Distance between targets (deg)');
ylabel('dprime (a.u.)');
title('Effect of distance when duration=0.5 s');
axis([0 0.75*pi 0 2]);
set(gca,'XTick',[0:pi/8:0.75*pi],'XTickLabel',(0:pi/8:0.75*pi)*180/pi);

drawPublishAxis('figSize=[10,5]');
savepdf(h,fullfile('~/proj/afcom/figures/distance.pdf'));
