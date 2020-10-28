function [fx,fxp,fy,fyp,...
            fypyp,fypy,fypxp,fypx,fyyp,fyy,fyxp,...
            fyx,fxpyp,fxpy,fxpxp,fxpx,fxyp,fxy,fxxp,fxx,...
            f] = FisherIST(isLog)

if nargin == 0
    isLog = 0;
end

syms par_beta par_delta par_gamma par_alpha par_nu
syms con k h invt out a v
syms conp kp hp invtp outp ap vp

f1 = exp((a + par_alpha*v)/(1-par_alpha))*con^(-1) - ....
    par_beta*(par_alpha*exp(ap)*kp^(par_alpha-1)*hp^(1-par_alpha) + ...
        (1-par_delta)*exp(-vp))*conp^(-1);
f2 = con - (1-par_alpha)*exp(a)*k^(par_alpha)*h^(-par_alpha);
f3 = con + invt - out;
f4 = out - exp(a)*k^(par_alpha)*h^(1-par_alpha);
f5 = exp((a+par_alpha*v)/(1-par_alpha)) * kp - exp(-v)*(1-par_delta)*k - invt;
f6 = ap - par_gamma;
f7 = vp - par_nu;

f = [f1; f2; f3; f4; f5; f6; f7];

x = [k a v];
y = [con h invt out];
xp = [kp ap vp];
yp = [conp hp invtp outp];

if isLog == 1
    f = subs(f, [k, con, invt, out, kp, conp, invtp, outp], exp([k, con, invt, out, kp, conp, invtp, outp]));
end
[fx,fxp,fy,fyp,...
    fypyp,fypy,fypxp,fypx,fyyp,fyy,fyxp,...
    fyx,fxpyp,fxpy,fxpxp,fxpx,fxyp,fxy,fxxp,fxx]=anal_deriv(f,x,y,xp,yp);
