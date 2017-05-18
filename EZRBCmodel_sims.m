EZRBCmodel_run;

sig = 0.01;
num_sim = 1000;
e = randn(num_sim, 1);


% isLog = 0;
% [fx,fxp,fy,fyp,fypyp,fypy,fypxp,fypx,fyyp,fyy,fyxp,fyx,fxpyp,fxpy,fxpxp,fxpx,fxyp,fxy,fxxp,fxx,f] = EZRBCmodel(isLog);
% [par_beta, par_delta, par_rho, par_gamma, par_alpha, par_lambda, par_phi, eta, con, invt, out, k, z, n, J, R, conp, invtp, outp, kp, zp, np, Jp, Rp] = EZRBCmodel_ss(isLog);
% approx = 2;
% num_eval;
% [gx,hx] = gx_hx(nfy,nfx,nfyp,nfxp);
% [gxx,hxx] = gxx_hxx(nfx,nfxp,nfy,nfyp,nfypyp,nfypy,nfypxp,nfypx,nfyyp,nfyy,nfyxp,nfyx,nfxpyp,nfxpy,nfxpxp,nfxpx,nfxyp,nfxy,nfxxp,nfxx,hx,gx);
% [gss,hss] = gss_hss(nfx,nfxp,nfy,nfyp,nfypyp,nfypy,nfypxp,nfypx,nfyyp,nfyy,nfyxp,nfyx,nfxpyp,nfxpy,nfxpxp,nfxpx,nfxyp,nfxy,nfxxp,nfxx,hx,gx,gxx,eta);

x0 = [k z];
[Y, X] = simu_2nd(gx, hx, gxx, hxx, gss, hss, eta, sig, x0, e);
figure(1); plot(X(:,1)); title('Capital') % capital
% figure(2); plot(X(:,2)); % exogenous shock
% figure(3); plot(Y(:,1)); % J
% figure(4); plot(Y(:,2)); % R
figure(5); plot(Y(:,3)); title('Consumption') % consumption
figure(6); plot(Y(:,4)); title('Labor') % labor
figure(7); plot(Y(:,5)); title('Investment') % investment
figure(8); plot(Y(:,6)); title('Output') % output

