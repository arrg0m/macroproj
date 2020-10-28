function [fx,fxp,fy,fyp,fypyp,fypy,fypxp,fypx,fyyp,fyy,fyxp,fyx,fxpyp,fxpy,fxpxp,fxpx,fxyp,fxy,fxxp,fxx,f] = Tallarini_model
syms 	SIG DELTA ALFA BETTA RHO
syms c cp k kp a ap 

f1 = c + kp - (1-DELTA) * k - a * k^ALFA;
f2 = c^(-SIG) - BETTA * cp^(-SIG) * (ap * ALFA * kp^(ALFA-1) + 1 - DELTA);
f3 = log(ap) - RHO * log(a);

f = [f1;f2;f3];

x = [k a];
y = c;
xp = [kp ap];
yp = cp;

f = subs(f, [x,y,xp,yp], (exp([x,y,xp,yp])));

[fx,fxp,fy,fyp,fypyp,fypy,fypxp,fypx,fyyp,fyy,fyxp,fyx,fxpyp,fxpy,fxpxp,fxpx,fxyp,fxy,fxxp,fxx]=anal_deriv(f,x,y,xp,yp);