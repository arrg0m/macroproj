function [par_beta, par_delta, par_gamma, par_alpha, par_nu, eta...
            con, k, h, invt, out, a, v, w,...
            conp, kp, hp, invtp, outp, ap, vp, wp] = FisherIST_ss(shock, isLog)
        
if nargin == 0
    shock = 2;
end
if nargin == 1
    isLog = 0;
end

par_beta = 0.99;
par_delta = 0.025;
par_gamma = 0.0026;
par_alpha = 0.33;
par_nu = 0.0046;
if shock == 0
    eta = [0 0; 1-par_alpha 0; 0 0];
elseif shock == 1
    eta = [0 0; 0 0; 0 1];
else
    eta = [0 0; 1-par_alpha 0; 0 1];
end

a = par_gamma;
v = par_nu;
khratio = ((par_alpha*par_beta*exp(par_gamma))/(exp((par_gamma+par_alpha*par_nu)/(1-par_alpha)) - par_beta*(1-par_delta)*exp(-par_nu)))^(1/(1-par_alpha));
con = (1-par_alpha)*exp(par_gamma)*khratio^par_alpha;
k = con/(exp(par_gamma)*khratio^(par_alpha-1) - exp((par_gamma+par_alpha*par_nu)/(1-par_alpha)) + (1-par_delta)*exp(-par_nu));
invt = (exp((par_gamma+par_alpha*par_nu)/(1-par_alpha)) - (1-par_delta)*exp(-par_nu))*k;
out = con + invt;
h = k/khratio;
w = (log(con) - h)/(1-par_beta);

if isLog == 1 % log-linearinzation
    con = log(con);
    k = log(k);
    invt = log(invt);
    out = log(out);
end

ap = a;
vp = v;
conp = con;
kp = k;
hp = h;
invtp = invt;
outp = out;
wp = w;