function [fx,fxp,fy,fyp,fypyp,fypy,fypxp,fypx,fyyp,fyy,fyxp,fyx,fxpyp,fxpy,fxpxp,fxpx,fxyp,fxy,fxxp,fxx,f] = EZGmodel(isLog)

if nargin == 0
    isLog = 0;
end

syms par_beta par_delta par_rho par_gamma par_alpha par_lambda par_phi
syms Js Rs Jn Rn cn cs z k n w invt out bs bn pf 
syms Jsp Rsp Jnp Rnp cnp csp zp kp np wp invtp outp bsp bnp pfp

f1 = Js^(1-par_rho) - (1-par_beta)*cs^(1-par_rho) - par_beta*Rs^(1-par_rho); 
f2 = Jn^(1-par_rho) - (1-par_beta)*(cn*(1-n)^(par_lambda))^(1-par_rho) - par_beta*Rn^(1-par_rho);
f3 = cs^(-par_rho) - par_beta*(par_alpha * zp * kp^(par_alpha-1)*np^(1-par_alpha) + 1 - par_delta) * csp^(-par_rho) * (Jsp/Rs)^(par_rho - par_gamma);
f4 = pf*cs^(-par_rho) - par_beta*csp^(-par_rho) * (Jsp/Rs)^(par_rho - par_gamma);
f5 = cs + invt + pf*bsp - out + w*n - bs;
f6 = out - z*k^(par_alpha)*n^(1-par_alpha);
f7 = kp - (1-par_delta)*k - invt;
f8 = Rs^(1-par_gamma) - Jsp^(1-par_gamma);
f9 = log(zp) - par_phi*log(z);
f10 = par_lambda*cn - w*(1-n);
f11 = pf*cn^(-par_rho) - par_beta*cnp^(-par_rho) * (Jnp/Rn)^(par_rho-par_gamma) * ((1-np)/(1-n))^(par_lambda*(1-par_rho));
f12 = cn + pf*bnp - w*n - bn;
f13 = Rn^(1-par_gamma) - Jnp^(1-par_gamma);
%f14 = cs + cn + invt - out;
f15 = bs + bn;
f16 = w*out - (1-par_alpha)*n;

f = [f1; f2; f3; f4; f5; f6; f7; f8; f9; f10; f11; f12; f13; f15; f16];

x = [k z];
y = [Js Rs Jn Rn cs cn n w invt out bs bn pf];
xp = [kp zp];
yp = [Jsp Rsp Jnp Rnp csp cnp np wp invtp outp bsp bnp pfp];

if isLog == 1
    f = subs(f, [x, cn, cs, n, w, invt, out, pf, xp, cnp, csp, np, wp, invtp, outp, pfp], ...
        (exp([x, cn, cs, n, w, invt, out, pf, xp, cnp, csp, np, wp, invtp, outp, pfp])));
end
[fx,fxp,fy,fyp,fypyp,fypy,fypxp,fypx,fyyp,fyy,fyxp,fyx,fxpyp,fxpy,fxpxp,fxpx,fxyp,fxy,fxxp,fxx]=anal_deriv(f,x,y,xp,yp);
