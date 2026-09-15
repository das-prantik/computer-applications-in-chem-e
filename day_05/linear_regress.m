function [slope, intercept, R2] = linear_regress(x, y) %#ok<FNDEF>

% @parameter x, y: data vectors

% @return slope, intercept: ordinary least-squares fit, y = intercept + slope*x
% @return R2: coefficient of determination

N = length(x);
Sx = sum(x); Sy = sum(y);
Sxx = sum(x .* x); Sxy = sum(x .* y);

slope = (N*Sxy - Sx*Sy) / (N*Sxx - Sx*Sx);
intercept = (Sy - slope*Sx) / N;

ybar = Sy / N;
ypred = intercept + slope*x;
SSt = sum((y - ybar).^2);
SSr = sum((y - ypred).^2);
R2 = 1 - SSr/SSt;

end
