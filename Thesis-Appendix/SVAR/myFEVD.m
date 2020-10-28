function Omega = myFEVD(Theta, is_diff)
    [k, kk, h] = size(Theta);
    if nargin == 2
        for itr = 1:k
            if is_diff(itr) == 1
                Theta(itr, :, :) = cumsum(Theta(itr, :, :), 3);
            end
        end
    end
    Thetasq = Theta.^2;
    MSEcomp = cumsum(Thetasq, 3);
    MSEsum = sum(MSEcomp, 2);
    Omega = MSEcomp./repmat(MSEsum, [1 k]);