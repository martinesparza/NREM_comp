function [dydx,k] = myode_two(x, y)
k = x.^2 + y.^2;
dydx = x + k.*y;

end

