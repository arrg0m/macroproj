function [par_beta, par_delta, par_gamma, par_alpha, eta...
            con, k, h, invt, out, a, ...
            conp, kp, hp, invtp, outp, ap] = FisherNshock_ss(isLog)
        
if nargin == 0
    isLog = 0;
end

par_beta = 0.99;
par_delta = 0.025;
par_gamma = 0.0026;
par_alpha = 0.33;
eta=[0 1]'; % error

a = par_gamma;
khratio = (par_beta*par_alpha*exp(par_gamma)/(exp(par_gamma/(1-par_alpha))-par_beta*(1-par_delta)))^(1/(1-par_alpha));
con = (1-par_alpha)*exp(par_gamma)*khratio^(par_alpha);
k = ((1-par_alpha)*exp(par_gamma)*khratio^(par_alpha))/(exp(par_gamma)*khratio^(par_alpha-1) - 1 + (1-par_delta));
h = k/khratio;
invt = (1-(1-par_delta)*exp(-par_gamma/(1-par_alpha)))*k;
out = exp(-par_gamma/(1-par_alpha))*con + invt;

if isLog == 1
    % log-linearinzation
    con = log(con);
    k = log(k);
    h = log(h);
    invt = log(invt);
    out = log(out);
end
   
ap = a;
conp = con;
kp = k;
hp = h;
invtp = invt;
outp = out;