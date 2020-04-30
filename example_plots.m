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
angles = -pi:pi/12:2*pi;
Lab = [repmat(L,length(angles),1) 128*cos(angles') 128*sin(angles')];
cform = makecform('lab2srgb');
rgb = applycform(Lab,cform);

rads = -pi:pi/64:pi;

h = figure; hold on

for ai = 1:length(angles)
    % shift rads by the angle
    rads_ = rads - angles(ai);
    rads_(rads_<0) = rads_(rads_<0)+2*pi;
    rads_(rads_>2*pi) = rads_(rads_>2*pi)-2*pi;
    % get the pscale values
    activation = 1 - pscale(abs(rads_*180/pi));
    % plot
    plot(rads,activation);
end

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
