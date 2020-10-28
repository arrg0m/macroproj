function [par_beta, par_delta, par_gamma, par_alpha, par_nu, eta...
            con, k, h, invt, out, a, v, ...
            conp, kp, hp, invtp, outp, ap, vp] = FisherIST_ss(isLog)
        
if nargin == 0
    isLog = 0;
end

par_beta = 0.99;
par_delta = 0.025;
par_gamma = 0.0026;
par_alpha = 0.33;
par_nu = 0.0046; % -0.0046
eta=[0 0; 1-par_alpha 0; 0 1]; % error

a = par_gamma;
v = par_nu;
khratio = ((par_alpha*par_beta*exp(par_gamma))/(exp((par_gamma+par_alpha*par_nu)/(1-par_alpha)) - par_beta*(1-par_delta)*exp(-par_nu)))^(1/(1-par_alpha));
con = (1-par_alpha)*exp(par_gamma)*khratio^par_alpha;
k = con/(exp(par_gamma)*khratio^(par_alpha-1) - exp((par_gamma+par_alpha*par_nu)/(1-par_alpha)) + (1-par_delta)*exp(-par_nu));
invt = (exp((par_gamma+par_alpha*par_nu)/(1-par_alpha)) - (1-par_delta)*exp(-par_nu))*k;
out = con + invt;
h = k/khratio;

if isLog == 1
    % log-linearinzation
    con = log(con);
    k = log(k);
    invt = log(invt);
    out = log(out);
end
% 
% a = 0.0026;
% v = -0.0046;
% con = 0.749015;
% k = 3.3317;
% h = -0.143731;
% invt = -0.479525;
% out = 1.00576;
% 
% if isLog == 0
%     con = exp(con);
%     k = exp(k);
%     h = exp(h);
%     invt = exp(invt);
%     out = exp(out);
% end

ap = a;
vp = v;
conp = con;
kp = k;
hp = h;
invtp = invt;
outp = out;