function [fx,fxp,fy,fyp,fypyp,fypy,fypxp,fypx,fyyp,fyy,fyxp,fyx,fxpyp,fxpy,fxpxp,fxpx,fxyp,fxy,fxxp,fxx,f] = FisherNshock(isLog)

if nargin == 0
    isLog = 0;
end

syms par_beta par_delta par_gamma par_alpha
syms con k h invt out a
syms conp kp hp invtp outp ap

f1 = exp(a/(1-par_alpha))*con^(-1) - par_beta*(par_alpha*exp(ap)*kp^(par_alpha-1)*hp^(1-par_alpha)+(1-par_delta))*conp^(-1);
f2 = con - (1-par_alpha)*exp(a)*k^(par_alpha)*h^(-par_alpha);
f3 = exp(-a/(1-par_alpha))*con + invt - out;
f4 = out - exp(-par_alpha*a/(1-par_alpha))*k^(par_alpha)*h^(1-par_alpha);
f5 = kp - (1-par_delta)*exp(-a/(1-par_alpha))*k - invt;
f6 = a - par_gamma;

f = [f1; f2; f3; f4; f5; f6];

x = [k a];
y = [con h invt out];
xp = [kp ap];
yp = [conp hp invtp outp];

if isLog == 1
    f = subs(f, [k, y, kp, yp], exp([k, y, kp, yp]));
end
[fx,fxp,fy,fyp,fypyp,fypy,fypxp,fypx,fyyp,fyy,fyxp,fyx,fxpyp,fxpy,fxpxp,fxpx,fxyp,fxy,fxxp,fxx]=anal_deriv(f,x,y,xp,yp);
