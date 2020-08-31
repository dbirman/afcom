function fit = ac_fitSatExpBeta( x, y)
%FITSATEXP Fit a saturating exponential to the power data
%%
a = 0.5;
k = 20;

base = 8;

[bestp] = lsqnonlin(@(p) satExpResid(x,y,p),[k base]);

x_ = x(1):.001:x(end);
y_ = 0.5 + 0.5-0.5*bestp(2).^(-x_*bestp(1));


fit.x = x;
fit.x_ = x_;
fit.y = y;
fit.y_ = y_;
fit.a = 0.5;
fit.k = bestp(1);
fit.base = bestp(2);
% fit.k = bestp(2);

function res = satExpResid(x,y,p)

k = p(1);
base = p(2);

res = y - (0.5 + 0.5-0.5*base.^(-x*k));