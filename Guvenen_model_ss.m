function [par_beta, par_delta, par_rho, par_gamma, par_alpha, par_lambda, par_phi,...
    pf, z, w, bs, bn, n, k, invt, out, cs, cn, ...
    pfp, zp, wp, bsp, bnp, np, kp, invtp, outp, csp, cnp] = Guvenen_model_ss

par_beta = 0.99;
par_delta = 0.025;
par_rho = 0.95;
par_gamma = 5;
par_alpha = 0.36;
par_lambda = 0.5; % 0.5, 1, or 2. TODO include as a funtional input
par_phi = 0.80; %???

pf = par_beta;
z = 1;
knratio = ((1/par_beta - (1-par_delta))/par_alpha)^(1/(par_alpha-1));
w = (1-par_alpha)*z*knratio^(par_alpha);
bs = 0.01; % IN FACT; CANNOT LOG-LINEARIZE at zero; or ...
bn = 0.01;
n = 1/(par_lambda + 1);
k = knratio * n;
invt = par_delta*k;
out = z*knratio*n;
cs = out - w*n - invt;
cn = w*n;

% pf = log(pf);
% z = log(z);
% w = log(w);
% bs = log(bs);
% bn = log(bn);
% n = log(n);
% k = log(k);
% invt = log(invt);
% out = log(out);
% cs = log(cs);
% cn = log(cn);

pfp = pf;
zp = z;
wp = w;
bsp = bs;
bnp = bs;
np = n;
kp = k;
invtp = invt;
outp = out;
csp = cs;
cnp = cs;