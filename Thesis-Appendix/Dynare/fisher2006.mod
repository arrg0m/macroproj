% U(C,H) = ln(C) - H

var y x a v c h k U W;
predetermined_variables k;
varexo ea ev;
parameters alpha delta beta nu gamma;

alpha = 0.33;
%delta = 1;
delta = 0.025;
beta = 0.99;
nu = -0.0046;
%nu = 0.0046;
gamma = 0.0026;

sigmaev = 0.01;
sigmaea = 0.01*(1-alpha);

model;
    exp(y) = exp(alpha*k + (1-alpha)*h + a);
    exp(c) + exp(x) = exp(y);
    exp((a+alpha*v)/(1-alpha) + k(+1)) = (1-delta)*exp(-v + k)+ exp(x);
    exp((a+alpha*v)/(1-alpha) - c) = beta*exp(-c(+1)) * (alpha*exp((alpha-1)*k(+1) + (1-alpha)*h(+1) + a(+1)) +(1-delta)*exp(-v(+1)));
    exp(c) = (1-alpha)*exp(alpha*k - alpha*h + a);
    U = c - exp(h);
    W = c - exp(h) + beta*W(+1);
    a = gamma + ea;
    v = nu + ev;
end;

initval;
    y = log(4);
    k = log(30);
    c = log(4);
    x = log(2);
    v = 0;
    a = 0;
    h = log(0.8);
    U = 0;
    W = 0;
end;
shocks;
    var ea = sigmaea^2;
    var ev = sigmaev^2;
end;
steady;

H = 40;
stoch_simul(hp_filter = 1600, irf = 40, order = 2);


% dynare fisher2006.mod