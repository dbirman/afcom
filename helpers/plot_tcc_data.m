function plot_tcc_data(xs,data,col,name)

if size(data,1)>1
    y = mean(data); % don't worry about error bars for now
else
    y = data;
end

%%
% if the xs are only >0, then flip them
if ~any(xs<0)
    xs = [fliplr(-xs) xs];
    y = [fliplr(y) y];
end

h = figure; hold on

r = 5;
off = pi/2;

circles(0,0,1,'facecolor','none','edgecolor',col);
for xi = 1:length(xs)
    % plot a vector in this direciton
    x1 = 1 * cos(xs(xi)+off);
    y1 = 1 * sin(xs(xi)+off);

    x2 = (1+r*y(xi))*cos(xs(xi)+off);
    y2 = (1+r*y(xi))*sin(xs(xi)+off);
    
    plot([x1 x2],[y1 y2],'-','LineWidth',1,'Color',col);
end

plot([0 0],[0 1],'--r');
axis([-1 1 -1 1]*3);
axis square

set(gca,'XTick',[0 1],'XTickLabel',{'',''});
set(gca,'YTick',[0 1],'YTickLabel',{'',''});


savepdf(h,fullfile('~/proj/afcom/figures/',name));