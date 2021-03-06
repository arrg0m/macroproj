function [par_beta, par_delta, par_gamma, par_alpha, par_phi,...
    pf, z, w, bs, bn, n, k, invt, out, cs, cn, Js, Rs, Jn, Rn,...
    pfp, zp, wp, bsp, bnp, np, kp, invtp, outp, csp, cnp, Jsp, Rsp, Jnp, Rnp] = EZGmodel_ss(isLog, par_lambda)

if nargin == 0
    isLog = 0;
end

par_beta = 0.99;
par_delta = 0.025;
par_gamma = 5;
par_alpha = 0.36;
par_phi = 0.95;


pf = par_beta;
z = 1;
knratio = ((1/par_beta - (1-par_delta))/par_alpha/z)^(1/(par_alpha-1));
w = (1-par_alpha)*z*knratio^(par_alpha);
n = 1/(par_lambda + 1);
k = knratio * n;
invt = par_delta*k;
out = z*knratio*n;
cn = w*(1-n)/par_lambda;
cs = out - cn - invt;
Js = cs;
Rs = Js;
Jn = cn*(1-n)^par_lambda;
Rn = Jn;
bn = (cn-w*n)/(1-pf);
bs = -bn;

if isLog == 1
    pf = log(pf);
    z = log(z);
    w = log(w);
    n = log(n);
    k = log(k);
    invt = log(invt);
    out = log(out);
    cs = log(cs);
    cn = log(cn);
end

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
Jsp = Js;
Rsp = Rs;
Jnp = Jn;
Rnp = Rn;