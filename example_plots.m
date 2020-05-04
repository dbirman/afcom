%% Compute colors for example plots

L = 60;
angles = 0:pi/128:2*pi;
Lab = [repmat(L,length(angles),1) 128*cos(angles') 128*sin(angles')];
cform = makecform('lab2srgb');
rgb = applycform(Lab,cform);

% Plot in a ring
h = figure; hold on
din = 3;
dout = 4;
for ai = 1:(length(angles)-1)
    astart = angles(ai);
    aend = angles(ai+1);
    
    % compute the inner and outer x positions
    xs = [din*cos(astart) din*cos(aend) dout*cos(aend) dout*cos(astart)];
    ys = [din*sin(astart) din*sin(aend) dout*sin(aend) dout*sin(astart)];
    
    p = fill(xs,ys,rgb(ai,:),'EdgeColor',rgb(ai,:));
%     set(p,'FaceColor','w');
end
axis square

drawPublishAxis('figSize=[4.5 4.5]');

savepdf(h,fullfile('~/proj/afcom/figures/colorwheel.pdf'));

%% Now using the same matched colors, generate the "channel" tunings

L = 60;
angles = -pi:pi/8:pi;
Lab = [repmat(L,length(angles),1) 128*cos(angles') 128*sin(angles')];
cform = makecform('lab2srgb');
rgb = applycform(Lab,cform);

rads = -pi:pi/128:pi;

h = figure; hold on

for ai = 1:length(angles)
    % shift rads by the angle
    rads_ = abs(rads - angles(ai));
%     rads_ = mod(rads_,pi);
    % get the pscale values
    activation = 1 - pscale(rads_*180/pi);
    activation = activation - min(activation);
    activation = activation ./ max(activation);
    % plot - set color according to rgb
    plot(rads,activation,'Color',rgb(ai,:));
end

axis([-pi pi 0 1]);
set(gca,'Xtick',[-pi 0 pi],'XTickLabel',[-180 0 180],'YTick',[0 1]);
xlabel('Angle (deg)');
ylabel('Channel activation (a.u.)');

drawPublishAxis('figSize=[3.5,3.5]');

savepdf(h,fullfile('~/proj/afcom/figures/colorchannels.pdf'));

%% Now plot the individual 
% channel responses for each of the four different dot patches, showing how
% the channels responded +/- 1 sigma noise
dps = [2 1 0.5 0.5];

targets = [3/4*pi 3/2*pi 1/4*pi 5/4*pi]-pi/2;

h = figure; hold on
for di = 1:length(dps)
    dp = dps(di);
    subplot(4,1,di); hold on
    % get the max activation for each channels
    dists = angdist(angles,targets(di));
    activations = dp * (1-pscale(abs(dists*180/pi)));
    % plot these
    for ai = 1:length(angles)
        errbar(angles(ai),activations(ai),1,'-','Color',rgb(ai,:));
        plot(angles(ai),activations(ai),'o','MarkerFaceColor',rgb(ai,:),'MarkerEdgeColor','w','MarkerSize',4);
%         plot(angles(ai),randn+activations(ai),'o','MarkerFaceColor','k','MarkerEdgeColor','w','MarkerSize',3);
    end

	

    axis([-pi pi -2 3.5]);
    
    set(gca,'Xtick',[-pi 0 pi],'XTickLabel',[-180 0 180],'YTick',[-2 -1 0 1 2]);
    
    drawPublishAxis('figSize=[3.5,8.9]');
end

savepdf(h,fullfile('~/proj/afcom/figures/coloractivations.pdf'));

%% Now plot the individual 
% channel responses for each of the four different dot patches, showing how
% the channels responded +/- 1 sigma noise
dps = [2 1 0.5 0.5];

targets = [3/4*pi 3/2*pi 1/4*pi 5/4*pi]-pi/2;

h = figure; hold on
for di = 1:length(dps)
    dp = dps(di);
    subplot(4,1,di); hold on
    % get the max activation for each channels
    dists = angdist(angles,targets(di));
    activations = dp * (1-pscale(abs(dists*180/pi)));
    % plot these
    for ai = 1:length(angles)
%         errbar(angles(ai),activations(ai),1,'-','Color',rgb(ai,:));
%         plot(angles(ai),activations(ai),'o','MarkerFaceColor',rgb(ai,:),'MarkerEdgeColor','w','MarkerSize',4);
        plot(angles(ai),randn+activations(ai),'o','MarkerFaceColor','k','MarkerEdgeColor','w','MarkerSize',3);
    end

	

    axis([-pi pi -2 3.5]);
    
    set(gca,'Xtick',[-pi 0 pi],'XTickLabel',[-180 0 180],'YTick',[-2 -1 0 1 2]);
    
    drawPublishAxis('figSize=[3.5,8.9]');
end

savepdf(h,fullfile('~/proj/afcom/figures/coloractivations_sampled.pdf'));

%% Now plot the TCC PDF function, rotated for each of the same target positions

h = figure;

angles = rads;
Lab = [repmat(L,length(angles),1) 128*cos(angles') 128*sin(angles')];
cform = makecform('lab2srgb');
rgb = applycform(Lab,cform);

likes = zeros(4,length(rads));

for ti = 1:length(targets)
    like = computeTCCPDF(rads,dps(ti));
    like(1) = like(2);
    like(end) = like(end-1);
    
    subplot(4,1,ti); hold on
    
    % circshift to target
    idx = find(rads>=targets(ti),1);
    like = circshift(like',idx+128);
    
    likes(ti,:) = like;
    
    for ri = 1:(length(rads)-1)
        rs = [rads(ri) rads(ri+1)];
        ys = [like(ri) like(ri+1)];
        plot(rs,ys,'-','Color',rgb(ri,:));
    end
    
    axis([-pi pi 0 0.06]);
    set(gca,'Xtick',[-pi 0 pi],'XTickLabel',[-180 0 180],'YTick',[0 0.02]);
    
    drawPublishAxis('figSize=[3,8.9]');
end

savepdf(h,fullfile('~/proj/afcom/figures/colorlikelihoods.pdf'));

%% And finally, plot the full likelihood
out = likes' * [0.5 0.25 0.2 0.05]';

h = figure; hold on

angles = rads;
Lab = [repmat(L,length(angles),1) 128*cos(angles') 128*sin(angles')];
cform = makecform('lab2srgb');
rgb = applycform(Lab,cform);

% instead of plotting the positions in black, let's color them to make it
% easier to see which color goes with which angle
for ri = 1:(length(rads)-1)
    rs = [rads(ri) rads(ri+1)];
    ys = [out(ri) out(ri+1)];
    plot(rs,ys,'-','Color',rgb(ri,:));
end
axis([-pi pi 0 0.025]);
set(gca,'Xtick',[-pi 0 pi],'XTickLabel',[-180 0 180],'YTick',[0 0.02]);

drawPublishAxis('figSize=[3,2]');

savepdf(h,fullfile('~/proj/afcom/figures/color_finallike.pdf'));
%% OLD CODE OLD CODE

%% example plot with different sensitivity (left panel scaled, right panel normalized)
x = -pi:pi/128:pi;

h = figure(2); clf; hold on

% cmap = colorblindmap/255;

cmap = brewermap(13,'Reds');

ax = pi*0.75;
ay = 0.4;
d = 0.5;

subplot(121); hold on
% compute the pscale function for four different angs
dprimes = logspace(-1,0.5,6);
dprimes = fliplr(dprimes);
legs = {};
for di = 1:length(dprimes)
    dprime = dprimes(di);
    
    activation = dprime*(1-pscale(abs(x)*180/pi));
    
    plot(x,activation,'-','Color',cmap(14-di,:));
    legs{end+1} = sprintf('d'' %1.2f',dprime);
end
legend(legs);


ylabel('Channel activation (sigma=1)');
xlabel('Distance from presented stimulus (deg)');
axis([-pi pi 0 3.2]);
set(gca,'XTick',[-pi 0 pi],'XTickLabel',[-180 0 180]);
set(gca,'YTick',[0 3.2]);
drawPublishAxis('figSize=[4.5,4.5]','poster=0');

subplot(122); hold on
% compute the pscale function for four different angs
dprimes = logspace(-1,0.5,6);
dprimes = fliplr(dprimes);
legs = {};
for di = 1:length(dprimes)
    dprime = dprimes(di);
    
    activation = computeTCCPDF(x,dprime);
    plot(x,activation,'-','Color',cmap(14-di,:));
    legs{end+1} = sprintf('d'' %1.2f',dprime);
end

ylabel('Response likelihood (pdf)');
xlabel('Distance from presented stimulus (deg)');
axis([-pi pi 0 0.08]);
set(gca,'XTick',[-pi 0 pi],'XTickLabel',[-180 0 180]);
set(gca,'YTick',[0 0.08]);
drawPublishAxis('figSize=[4.5,4.5]','poster=0');

savepdf(h,fullfile('~/proj/afcom/figures/','sensitivity_example.pdf'));

%%
x = -pi:pi/64:pi;
% also plot for these dprime values the tcc_data plots
for di = 1:length(dprimes)
    like = computeTCCPDF(x,dprimes(di));
    plot_tcc_data(x,like*4,cmap(14-di,:),sprintf('tcc_circle_%i.pdf',di));
end
