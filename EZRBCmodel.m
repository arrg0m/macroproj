function [fx,fxp,fy,fyp,fypyp,fypy,fypxp,fypx,fyyp,fyy,fyxp,fyx,fxpyp,fxpy,fxpxp,fxpx,fxyp,fxy,fxxp,fxx,f] = EZRBCmodel(isLog)


if nargin == 0
    isLog = 0;
end


syms par_beta par_delta par_rho par_gamma par_alpha par_lambda par_phi
syms k z J R con n invt out
syms kp zp Jp Rp conp np invtp outp

f1 = J^(1-par_rho) - (1-par_beta)*(con*(1-n)^(par_lambda))^(1-par_rho) + par_beta*R^(1-par_rho);
f2 = R^(1-par_gamma) - Jp^(1-par_gamma);
f3 = par_lambda*con - (1-par_alpha)*z*k^(par_alpha)*n^(-par_alpha)*(1-n);
f4 = (1-n)^(par_lambda*(1-par_rho))*con^(-par_rho) - par_beta*(par_alpha*zp*kp^(par_alpha-1)*np^(1-par_alpha) + 1 - par_delta)*(1-np)^(par_lambda*(1-par_rho))*conp^(-par_rho)*(Jp/R)^(par_rho-par_gamma);
f5 = con + invt - out;
f6 = out - z*k^(par_alpha)*n^(1-par_alpha);
f7 = kp - (1-par_delta)*k - invt;
f8 = log(zp) - par_phi*log(z);

f = [f1; f2; f3; f4; f5; f6; f7; f8];

x = [k z];
y = [J R con n invt out];
xp = [kp zp];
yp = [Jp Rp conp np invtp outp];

if isLog == 1
    f = subs(f, [x, con, n, invt, out, xp, conp, np, invt, outp], (exp([x, con, n, invt, out, xp, conp, np, invt, outp])));
end
[fx,fxp,fy,fyp,fypyp,fypy,fypxp,fypx,fyyp,fyy,fyxp,fyx,fxpyp,fxpy,fxpxp,fxpx,fxyp,fxy,fxxp,fxx]=anal_deriv(f,x,y,xp,yp);