function ac_plotSepLikes_tcc(fit)
folder = fit.dataType;
params = fit.params;
cmap = ac_cmap;

to = -0.025;
sgap = 0.005;
bgap = 0.01;
%% go through the five distributions and plot them on a figure
% scale the distributions appropriately
maxY = 0.1;
groups = {'cue_4','cue_side','cue_feat','cue_target','cue_1'};
for tt = 1:5
    
    h = figure;

    cmap = ac_cmap;
    targets = {'target','side','feat','dist'};
    targetNames = {'Target','Same side','Same feature','Distractor'};
    for ti = 1:length(targets)
        % plot target
        subplot(3,4,ti); hold on

        if tt<5 && ti==1
            % also plot the cue_1 condition on top of this
            plot(fit.x,squeeze(fit.out(5,ti,:)),'-','Color',[0.75 0.75 0.75]);
        end
        plot(fit.x,squeeze(fit.out(tt,ti,:)),'-','Color',cmap.(targets{ti}));
        
        
        axis([-pi pi 0 maxY]);
        set(gca,'XTick',[-pi 0 pi],'XTickLabel',[-180 0 180]);
        set(gca,'YTick',[0 maxY],'YTickLabel',{'',''});
        title(targetNames{ti});
        if ti==1
            ylabel('Response likelihood (pdf)');
            xlabel('Distance from stimulus angle');
        end
        drawPublishAxis('poster=1');
    end
    
    bs = fit.params.(sprintf('bs_%i',tt));
    bf = fit.params.(sprintf('bf_%i',tt));
    bi = fit.params.(sprintf('bi_%i',tt));
    betas = [bs*bf bs*(1-bf) (1-bs)*(1-bi) (1-bs)*bi];
    
    subplot(3,4,5:8); hold on
    % plot the proportion that are on the correct side
    rectangle('Position',[0 0 bs-bgap 1],'FaceColor',[0.75 0.75 0.75],'EdgeColor','w');
    plot([0 0 bs-bgap bs-bgap],[0 1 1 0],'-k');
    text(to+bs/2,0.5,sprintf('%2.0f%%',bs*100));
    text(to+bs/2,1.5,'Correct side');
    
    plot([bs+bgap bs+bgap 1 1],[0 1 1 0],'-k');
    text(to+1-(1-bs)/2,0.5,sprintf('%2.0f%%',(1-bs)*100));
    text(to+1-(1-bs)/2,1.5,'Wrong');
    
    axis([0 1 -1 2]);
    set(gca,'XTick',[0 1],'XTickLabels',{'',''});
    set(gca,'YTick',[0 1],'YTickLabels',{'',''});
    
    drawPublishAxis('poster=1');
    
    
    subplot(3,4,9:12); hold on
    
    % plot the side: (target and side)
    target = bs*bf;
    side = bs*(1-bf);
    rectangle('Position',[0 0 target-sgap 1],'FaceColor',cmap.target,'EdgeColor','w');
    plot([0 0 target-sgap target-sgap],[0 1 1 0],'-k');
    text(to+target/2,0.5,sprintf('%2.0f%%',target*100));
    text(to+target/2,1.5,'Target');
    
    if side>0.01
        rectangle('Position',[target+sgap 0 side-bgap-sgap 1],'FaceColor',cmap.side,'EdgeColor','w');
        plot([target+sgap target+sgap target+side-bgap target+side-bgap],[0 1 1 0],'-k');
        if side >0.04
            text(to+target+side/2,0.5,sprintf('%2.0f%%',side*100));
            text(to+target+side/2,1.5,'Same side');
        end
    end
    
    % plot the offsides: (feature and distractor)
    feat = (1-bs)*(1-bi);
    dist = (1-bs)*bi;
    if feat>0.01
        rectangle('Position',[target+side+bgap 0 feat-sgap-bgap 1],'FaceColor',cmap.feat,'Edgecolor','w');
        plot([target+side+bgap target+side+bgap target+side+feat-sgap target+side+feat-sgap],[0 1 1 0],'-k');
        if feat >0.04
            text(to+target+side+feat/2,0.5,sprintf('%2.0f%%',feat*100));
            text(to+target+side+feat/2,1.5,'Same feature');
        end
    end
%     
    if dist>0.01
        rectangle('Position',[target+side+feat+sgap 0 1-(target+side+feat+sgap) 1],'FaceColor',cmap.dist,'Edgecolor','w');
        plot([target+side+feat+sgap target+side+feat+sgap 1 1],[0 1 1 0],'-k');
    %     text(to+1-(1-side)/2,0.5,sprintf('%2.0f%%',(1-side)*100));
    %     text(to+1-(1-side)/2,1.5,'1-\beta_{side}');
    end
    
    axis([0 1 -1 2]);
    set(gca,'XTick',[0 1],'XTickLabel',{'',''});
    set(gca,'YTick',[0 1],'YTickLabel',{'',''});
    drawPublishAxis('figSize=[30,20]','poster=1');
    
    savepdf(h,fullfile('~/proj/afcom/figures',folder,sprintf('tcc_%s.pdf',groups{tt})));
end
