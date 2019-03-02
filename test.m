x = 0:5:360;
d = pscale(abs(x-180));
% d = d./max(d);

dspace = 0:.1:3;

like = zeros(length(x),length(dspace));

for xi = 1:length(x)
    for di = 1:length(dspace)
        like(xi,di) = normpdf(dspace(di),1-d(xi),1);
    end
end

h = figure(1);
imagesc(like);

% now compute the true likelihood function:
%  which is the probability that a given PDF is greater than ALL OTHER PDFs
%  we will do this in log likelihood space... 
dprimes = logspace(0,1,10);
% dprimes = 10;
like_ = zeros(length(x),length(dprimes));
legends = {};
for di = 1:length(dprimes)
    dprime = dprimes(di);
    for xi = 1:length(x)
        dh = 1-d(xi);
        for ci = 1:length(x)
            if xi~=ci
                dl = 1-d(ci);
                % not the same row
                like_(xi,di) = like_(xi,di) + log(normcdf(dprime*(dh-dl),0,1));
                % equivalent to
%                 like_(xi,di) = like_(xi,di) + log(normcdf(dh-dl,0,dprime));
            end
        end
    end
    legends{end+1} = sprintf('d'' = %1.2f',dprime);
end

%%
h = figure(2);
clf
hold on
plot(like_);
legend(legends);

%% add a single gaussian
figure;
hold on
y = like_(:,1);
y = y-min(y);
y = y/sum(y);
plot(x,y);

y2 = 1-pscale(abs(x-180));
y2 = y2./sum(y2);
plot(x,y2,'-r');

%% normalize them, to see if they're the same
figure;
hold on
for di = 1:length(dprimes)
    temp = like_(:,di);
    temp = temp./max(temp);
    plot(temp);
end

%%
y1 = 1-pscale(abs(x-180)); 

y2 = y1.^0.1;%normpdf(x,180,60);
% y2 = y2./max(y2);

figure; hold on
plot(x,y1,'-b');
% plot(x,y2,'-r');