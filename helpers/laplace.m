function y = laplace(x,mu,b)

y = zeros(size(x));

y(x<mu) = 1/2/b * exp(-(mu-x(x<mu))/b);
y(x>=mu) = 1/2/b * exp(-(x(x>=mu)-mu)/b);