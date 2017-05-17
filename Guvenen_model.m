function [fx,fxp,fy,fyp,fypyp,fypy,fypxp,fypx,fyyp,fyy,fyxp,fyx,fxpyp,fxpy,fxpxp,fxpx,fxyp,fxy,fxxp,fxx,f] = Guvenen_model

syms par_beta par_delta par_rho par_gamma par_alpha par_lambda par_phi
syms cs cn bs bn k out z invt pf w n
syms csp cnp bsp bnp kp outp zp invtp pfp wp np

fs1 = cs^(-par_gamma) - par_beta*csp^(-par_gamma)*(par_alpha * zp * kp^(par_alpha-1) * np^(1-par_alpha) + (1-par_delta));
fs2 = pf*cs^(-par_gamma) - par_beta*csp^(-par_gamma);
fs3 = cs + invt + pf*bsp - out + w*n - bs;
fs4 = out - z * k^(par_alpha) * n^(1-par_alpha);
fs5 = kp - (1-par_delta)*k - invt;
fs6 = log(zp) - par_phi*log(z);
fw1 = par_lambda*cn - w*(1-n);
fw2 = pf*(1-n)^(par_lambda*(1-par_gamma)) * cn^(-par_gamma) - par_beta*(1-np)^(par_lambda*(1-par_gamma))*cnp^(-par_gamma);
fw3 = cn + pf*bnp - w*n + bn; % bn <= -bn
% fm1 = cs + cn + invt - out; % Note that this condition is redundant
fm2 = bs - bn; % bn <= -bn
fm3 = w*n - (1-par_alpha)*out;

f = [fs1; fs2; fs3; fs4; fs5; fs6; fw1; fw2; fw3; fm2; fm3];

x = [k z];
y = [cs cn bs bn out invt pf w n];
xp = [kp zp];
yp = [csp cnp bsp bnp outp invtp pfp wp np];

% f = subs(f, [x, y, xp, yp], (exp([x, y, xp, yp])));

[fx,fxp,fy,fyp,fypyp,fypy,fypxp,fypx,fyyp,fyy,fyxp,fyx,fxpyp,fxpy,fxpxp,fxpx,fxyp,fxy,fxxp,fxx]=anal_deriv(f,x,y,xp,yp);
