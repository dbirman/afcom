%%
dprimes = 0.25:0.5:1.5;

for di = 1:length(dprimes)
    dprime = dprimes(di);
    
    x = -90:2.5:90;
    fx = 1-pscale(abs(x));%normpdf(x,0,1);
    fx = fx * dprime;
    sigma = 0.5;%sqrt(1000);
    % plot(fx)

    range = linspace(min(fx)-4*sigma,max(fx)+4*sigma,30) ;

    px = nan(size(x));
    for idx1 = 1:length(x)
        px_idx1 = 0;
        for a = range
            cdf_all = normpdf(a,fx(idx1),sigma)*(range(2)-range(1));
            for idx2 = 1:length(x)
                if idx2 == idx1
                    continue
                else
                    cdf_all = cdf_all*normcdf(a,fx(idx2),sigma);
                end
            end
            px_idx1 = px_idx1+cdf_all;
        end
        px(idx1) = px_idx1;
    end
    
    ps{di} = px;
end

figure; hold on
for di = 1:length(ps)
    plot(ps{di});
end

%%
figure(1);
clf
hold on
plot(x,px);
% plot(x,cumsum(px));
y2 = (1-pscale(abs(x/sigma)));
y2 = y2 ./ sum(y2);
plot(x,y2,'-g');%cumsum(y2)/sum(y2));
% figure; hold on;plot((fx-min(fx))/sum(fx-min(fx)),'b');plot((px-min(px))/sum(px-min(px)),'r');
% figure; hold on;plot(fx/sum(fx),'b');plot(px,'r');
% 
% figure;subplot(2,1,1);plot(fx,'b');subplot(2,1,2);plot(px,'r');

%% 