folders = {'s030020181023'};

%% Save data (only run this if you have the raw data folders)
for ci = 1:length(folders)
    savedata_offset(folders{ci},1);
end

%% Run analysis

fi = 1;
load(fullfile('~/Box Sync/AFCOM_DATA/',sprintf('data_%s.mat',folders{fi})));

%% Plot pRF x coordinates in each condition

clear r2 ang ecc sz
for scan = 1:3
    r2(:,:,:,scan) = data.pRF.overlays(1).data{scan};
    ang(:,:,:,scan) = data.pRF.overlays(2).data{scan};
    ecc(:,:,:,scan) = data.pRF.overlays(3).data{scan};
    sz(:,:,:,scan) = data.pRF.overlays(4).data{scan};
end

%% convert to x/y
x = ecc .* cos(ang);
y = ecc .* sin(ang);

%% plot all x/y values from relevant indices
idx = r2>0.3;
idx = idx .* (sz<25);
idx = idx .* (x<15);
idx = idx .* (x>-15);
idx = idx .* (y<10);
idx = idx .* (y>-10);

keep = all(idx,4); % make sure all scans were relevant
disp(sum(keep(:)));

%% check quantities
clear cx cy
for scan = 1:3
    tx = x(:,:,:,scan);
    ty = y(:,:,:,scan);
    cx(:,scan) = tx(keep);
    cy(:,scan) = ty(keep);
end

%%
dLx = cx(:,2)-cx(:,1);
dLy = cy(:,2)-cy(:,1);
mu=mean(dLx);
ci=bootci(1000,@mean,dLx);
disp(sprintf('Left attend X shift: %1.2f, 95%% CI [%1.2f %1.2f]',mu,ci(1),ci(2)));
mu=mean(dLy);
ci = bootci(1000,@mean,dLy);
disp(sprintf('Left attend Y shift: %1.2f, 95%% CI [%1.2f %1.2f]',mu,ci(1),ci(2)));

dRx = cx(:,3)-cx(:,1);
dRy = cy(:,3)-cy(:,1);
mu = mean(dRx);
ci = bootci(1000,@mean,dRx);
disp(sprintf('Right attend X shift: %1.2f, 95%% CI [%1.2f %1.2f]',mu,ci(1),ci(2)));
mu = mean(dRy);
ci = bootci(1000,@mean,dRy);
disp(sprintf('Right attend Y shift: %1.2f, 95%% CI [%1.2f %1.2f]',mu,ci(1),ci(2)));

%% get only vertices 
figure; hold on
cmap = brewermap(3,'Dark2');
for c1 = 1:size(x,1)
    for c2 = 1:size(x,2)
        for c3 = 1:size(x,3)
            if keep(c1,c2,c3)
                cx = squeeze(x(c1,c2,c3,:));
                cy = squeeze(y(c1,c2,c3,:));

                % plot two lines, one from 1->2 and 1->3
                p1 = plot(cx(1:2),cy(1:2),'-','Color',cmap(2,:));
                p2 = plot(cx([1 3]),cy([1 3]),'-','Color',cmap(3,:));
            end
        end
    end
end
legend({'Attend left','Attend right'});

%% Plot arrows for V1
for roi = 8
    figure; hold on
    coords = data.coords{roi};
    
    for ci = 1:length(coords)
        coord = coords(:,ci);
        c1 = coord(1); c2 = coord(2); c3 = coord(3);
        if keep(c1,c2,c3)
            cx = squeeze(x(c1,c2,c3,:));
            cy = squeeze(y(c1,c2,c3,:));

            % plot two lines, one from 1->2 and 1->3
            arrow([cx(1) cy(1)],[cx(2) cy(2)],'Color',cmap(2,:));
            arrow([cx(1) cy(1)],[cx(3) cy(3)],'Color',cmap(3,:));
%             p1 = plot(cx(1:2),cy(1:2),'-','Color',cmap(2,:));
%             p2 = plot(cx([1 3]),cy([1 3]),'-','Color',cmap(3,:));
        end
    end
end