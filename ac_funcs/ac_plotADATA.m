function h = ac_plotADATA(~,adata,fit)
% Plot the adata structure showing the different conditions as histograms.
% Optional argument takes in a fit (from basic or encoding models) and
% plots those on top of the histograms.

% 
% headers = 
% 
%   Columns 1 through 7
% 
%     'target'    'trialType'    'cue'    'duration'    'dead'    'targetAngle'    'distractorAngle'
% 
%   Columns 8 through 13
% 
%     'angle1'    'angle2'    'angle3'    'angle4'    'respAngle'    'respDistance'
% 
%   Columns 14 through 16
% 
%     'distDistance'    'runs'    'RT'


cmap_ = colorblindmap/255;
global cmap
cmap = [];
cmap(5,:) = [0.5 0.5 0.5];
cmap(4,:) = cmap_(7,:);
cmap(3,:) = cmap_(4,:);
cmap(2,:) = cmap_(8,:);
cmap(1,:) = cmap_(1,:);

%% remove dead trials
adata = sel(adata,5,0);

%% split by cue
colordata = sel(adata,3,1);
dirdata = sel(adata,3,2);

%% for each trialType 0:4, pull the respAngle and targetAngle and plot

h = figure;


bins = -pi:pi/16:pi;


if ~isempty(fit)
    x = fit{1,1}.x;
end

durs = [1 0.25];

for tt = 0:4
    for di = 1:length(durs)
        dat = sel(colordata,2,tt);
        dat = sel(dat,4,durs(di));

        subplot(4,5,(di-1)*5+tt+1); hold on
        phelper_(dat,tt,bins);
        
        colorout = fit{1,di}.out;
        dirout = fit{2,di}.out;

        if ~isempty(fit)
            outx = colorout(tt+1,:);
            outb = interp1(x,outx,bins);
            outb = outb ./ sum(outb);
            plot(bins,outb,'-','Color',cmap(tt+1,:));
        end

        % SAME THING FOR DIRECTION
        dat = sel(dirdata,2,tt);
        dat = sel(dat,4,durs(di));
        
        subplot(4,5,10+(di-1)*5+tt+1); hold on

        phelper_(dat,tt,bins);

        if ~isempty(fit)
            outx = dirout(tt+1,:);
            outb = interp1(x,outx,bins);
            outb = outb ./ sum(outb);
            plot(bins,outb,'-','Color',cmap(tt+1,:));
        end
    end
end

function phelper_(dat,tt,bins)
%%
global cmap
hold on

groups = {'all','spatial','feature','target','baseline'};
% pull the respAngle and targetAngle
ta = dat(:,6);
ra = dat(:,12);

% rotate everything to a target angle of zero
ra = mod(ra-ta+pi,2*pi)-pi;

% to generate error bars we will take bootstraps over the dataset and
% re-compute the histogram multiples time
for i = 1:100 % number of bootstraps
    ra_boot = randsample(ra,length(ra),true);
    n(i,:) = hist(ra_boot,bins)./length(ra);
end

mu = mean(n);
ci = bootci(100,@mean,n);

errbar(bins,mu,ci(2,:)-mu,'Color',cmap(tt+1,:));
plot(bins,mu,'o','MarkerFaceColor',cmap(tt+1,:),'MarkerEdgeColor','w');

axis([-pi pi 0 0.4]);
set(gca,'XTick',[-pi -pi/2 0 pi/2 pi],'XTickLabel',{'-pi','-pi/2','0','pi/2','pi'});
xlabel('Response distance (degs)');
if tt==0
    ylabel('Proportion of trials (%)');
end

drawPublishAxis('figSize=[30,30]');