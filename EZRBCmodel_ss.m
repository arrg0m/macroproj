function [par_beta, par_delta, par_rho, par_gamma, par_alpha, par_lambda, par_phi, eta,...
            con, invt, out, k, z, n, J, R, ...
            conp, invtp, outp, kp, zp, np, Jp, Rp] = EZRBCmodel_ss(isLog)
        
if nargin == 0
    isLog = 0;
end

par_beta = 0.99;
par_delta = 0.025;
par_rho = 0.95;
par_gamma = 5;
par_alpha = 0.36;
par_lambda = 0.5; % 0.5, 1, or 2. TODO include as a funtional input
par_phi = 0.8; %???
eta=[0 1]';

z = 1;
knratio = ((1/par_beta - (1-par_delta))/par_alpha/z)^(1/(par_alpha - 1));
n = 1/(par_lambda*(1 + 1/par_lambda - knratio^(par_alpha-1)*par_delta/z ));
con = z*knratio^(par_alpha)*(1-n)/par_lambda;
out = z*knratio^(par_alpha)*n;
R = con*(1-n)^par_lambda;
J = R;
k = knratio*n;
invt = par_delta*k;

if isLog == 1
    % log-linearinzation
    z = log(z);
    n = log(n);
    con = log(con);
    out = log(out);
    k = log(k);
    invt = log(invt);
end
   
zp = z;
np = n;
conp = con;
outp = out;
Rp = R;
Jp = J;
kp = k;
invtp = invt;