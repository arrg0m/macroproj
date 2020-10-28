% HP Filter 
% Based on the formulation of http://www.auburn.edu/~hzk0001/hpfilter.pdf
function[trend, cycle] = HPfilter(y, lambda)
    y = y(:);
    len = size(y, 1);
    if len <= 5
        trend = y;
        cycle = 0*y;
    else
        F = diag([1; 5; 6*ones(len-4, 1); 5; 1]) + ...
            diag([-2; -4*ones(len-3, 1); -2], 1) + diag([-2; -4*ones(len-3, 1); -2], -1) + ...
            diag(ones(len-2, 1), 2) + diag(ones(len-2, 1), -2);
        trend = (lambda*F + eye(len))\y;
        cycle = y - trend;
    end