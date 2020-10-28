% IV Estimation
% 
% econometrics toolbox fgls ;(
% my fgls?
function[coef, resid] = myIV(y, X, Z, lambda)
    ZX = Z*(Z'*Z)^(-1)*Z'*X;
    Zy = Z*(Z'*Z)^(-1)*Z'*y;
    if nargin == 4
        coef = (ZX'*ZX + lambda*eye(size(ZX,2)))\ZX'*Zy;
    end
    resid = y - X*coef;
