function [fx,fxp,fy,fyp,fypyp,fypy,fypxp,fypx,fyyp,fyy,fyxp,fyx,fxpyp,fxpy,fxpxp,fxpx,fxyp,fxy,fxxp,fxx,f] = Nshock(isLog)


if nargin == 0
    isLog = 0;
end

syms par_beta par_delta par_rho par_gamma par_alpha par_lambda par_phi
syms 
syms 

f1 = con + invt - a * k^(par_alpha) * hr^(1-par_alpha);
f2 = kp - (1-par_delta)*k - invt;
f3 = 

f = [f1; f2; f3; f4; f5; f6; f7; f8];

x = [k z];
y = [J R con n invt out];
xp = [kp zp];
yp = [Jp Rp conp np invtp outp];

if isLog == 1
    f = subs(f, [x, con, n, invt, out, xp, conp, np, invt, outp], (exp([x, con, n, invt, out, xp, conp, np, invt, outp])));
end
[fx,fxp,fy,fyp,fypyp,fypy,fypxp,fypx,fyyp,fyy,fyxp,fyx,fxpyp,fxpy,fxpxp,fxpx,fxyp,fxy,fxxp,fxx]=anal_deriv(f,x,y,xp,yp);