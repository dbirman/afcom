clear likee xs
encs = 5:5:80;
for ei = 1:length(encs)
    enc = encs(ei);
    x = linspace(0,180,enc);
    d = 1-pscale(abs(x));
    
    % precompute normcdf
    % xn = -10:.0001:10;
    % y = normcdf(xn);

    rangeNeeded = [];

    % now compute the true likelihood function:
    %  which is the probability that a given PDF is greater than ALL OTHER PDFs
    %  we will do this in log likelihood space... 
    dprimes = 0.05; %logspace(0,0.5,10);
    % dprimes = 10;
    like_ = zeros(length(x),length(dprimes));
    legends = {};
    for di = 1:length(dprimes)
        dprime = dprimes(di);
        for xi = 1:length(x)

            like_(xi,di) = prod(normcdf(dprime*(d(xi)-d),0,1));
    %         
    %         like_(xi,di) = sum(log((normcdf(dprime*(d(xi)-d),0,1))));
        end
        legends{end+1} = sprintf('d'' = %1.2f',dprime);
    end
    xs{ei} = x;
    likee{ei} = like_;
end

%%
h = figure(2);
clf
hold on
% normalize
for e = 1:length(xs)
    x = xs{e};
    like_ = likee{e};
    like_ = like_ ./ repmat(sum(like_),size(like_,1),1);
    plot(x,like_);
    
end
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
    % move into positive space
    temp = temp - min(temp);
    temp = temp./max(temp);
    plot(temp);
end
legend(legends);

%%
y1 = 1-pscale(abs(x-180)); 

y2 = y1.^0.1;%normpdf(x,180,60);
% y2 = y2./max(y2);

figure; hold on
plot(x,y1,'-b');
% plot(x,y2,'-r');


%% Conditional probability version
% what is the probability that a draw a is greater than all the other
% draws, computed for every possible value of a

x = 0:5:360;
as = 0:.5:20;

d = 1-pscale(abs(x-180));

dprimes = logspace(0,1,10);
clear like
for di = 1:length(dprimes)
    dprime = dprimes(di);
    
    ap = normpdf(as,dprime,1);
    ap = ap ./ sum(ap);
    
    for xi = 1:length(x)
        nxi = setdiff(1:length(x),xi);
        
        likea = 0;
        for ai = 1:length(as)
            a = as(ai);
                            % P(X=A)
            likea = likea + ap(ai) * prod(normcdf(ap(ai),dprime*d(nxi),1));
        end
        
        like(xi,di) = likea;
    end
end

%%
figure;
plot(x,log(like));
