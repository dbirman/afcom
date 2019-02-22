function ac_plotSepLikes(fit)
folder = fit.dataType;
params = fit.params;

%% go through the five distributions and plot them on a figure
% scale the distributions appropriately

groups = {'cue_4','cue_side','cue_feat','cue_target','cue_1'};
cmap = ac_cmap;
for gi = 1:length(groups)
    kappa = params.(sprintf('kappa%i',gi));
%     ac_plotSingleKappa(kappa);
    lapse = params.(sprintf('lapse%i',gi));
    bs = params.(sprintf('beta_side%i',gi));
    bf = params.(sprintf('beta_feat%i',gi));
    bd = params.(sprintf('beta_dist%i',gi));
    
    x = -pi:pi/128:pi;
    y = vonMises(x,0,kappa);
    ysum = sum(y);
    y = y./ ysum;
    
%     maxY = max(y);
    maxY = 0.03;
    
    h = figure;
    
    % plot target
    subplot(1,5,5); hold on
    plot(x,y*(1-lapse)*bs*bf,'Color',cmap.target);
    axis([-pi pi 0 maxY]);
    set(gca,'XTick',[-pi 0 pi],'XTickLabel',[-180 0 180]);
    set(gca,'YTick',[0 maxY],'YTickLabel',{'',''});
    title('Target');
    ylabel('Response likelihood (pdf)');
    xlabel('Distance from stimulus angle');
    drawPublishAxis;
    % plot same-side distractor
    subplot(1,5,4);
    plot(x,y*(1-lapse)*bs*(1-bf),'Color',cmap.side);
    axis([-pi pi 0 maxY]);
    set(gca,'XTick',[-pi 0 pi],'XTickLabel',[-180 0 180]);
    set(gca,'YTick',[0 maxY],'YTickLabel',{'',''});
    title('Same-side');
    drawPublishAxis;
    % plot feature distractor
    subplot(1,5,3);
    plot(x,y*(1-lapse)*(1-bs)*bd,'Color',cmap.feat);
    axis([-pi pi 0 maxY]);
    set(gca,'XTick',[-pi 0 pi],'XTickLabel',[-180 0 180]);
    set(gca,'YTick',[0 maxY],'YTickLabel',{'',''});
    title('Same-feature');
    drawPublishAxis;
    % plot distractor distractor
    subplot(1,5,2);
    plot(x,y*(1-lapse)*(1-bs)*(1-bd),'Color',cmap.dist);
    axis([-pi pi 0 maxY]);
    set(gca,'XTick',[-pi 0 pi],'XTickLabel',[-180 0 180]);
    set(gca,'YTick',[0 maxY],'YTickLabel',{'',''});
    title('Irrelevant');
    drawPublishAxis;
    % plot lapse rate
    subplot(1,5,1);
    plot(x,repmat(vonMises(0,0,0)*lapse/ysum,1,length(x)),'Color',cmap.lapse);
    axis([-pi pi 0 maxY]);
    set(gca,'XTick',[-pi 0 pi],'XTickLabel',[-180 0 180]);
    set(gca,'YTick',[0 maxY],'YTickLabel',{'',''});
    title('Lapse');
    drawPublishAxis('figSize=[20 6]');
    
    savepdf(h,fullfile('~/proj/afcom/figures',folder,sprintf('allEncoding_%s.pdf',groups{gi})));
end