isLog = 1;
par_lambda = 0.5; % 0.5, 1, 2
par_rho = 1/0.05; % 1/1.5, 1/0.05;

EZRBCmodel_run;

sig = 0.01;
num_sim = 10000;
e = randn(num_sim, 1);

% isLog = 0;
% [fx,fxp,fy,fyp,fypyp,fypy,fypxp,fypx,fyyp,fyy,fyxp,fyx,fxpyp,fxpy,fxpxp,fxpx,fxyp,fxy,fxxp,fxx,f] = EZRBCmodel(isLog);
% [par_beta, par_delta, par_rho, par_gamma, par_alpha, par_lambda, par_phi, eta, con, invt, out, k, z, n, J, R, conp, invtp, outp, kp, zp, np, Jp, Rp] = EZRBCmodel_ss(isLog);
% approx = 2;
% num_eval;
% [gx,hx] = gx_hx(nfy,nfx,nfyp,nfxp);
% [gxx,hxx] = gxx_hxx(nfx,nfxp,nfy,nfyp,nfypyp,nfypy,nfypxp,nfypx,nfyyp,nfyy,nfyxp,nfyx,nfxpyp,nfxpy,nfxpxp,nfxpx,nfxyp,nfxy,nfxxp,nfxx,hx,gx);
% [gss,hss] = gss_hss(nfx,nfxp,nfy,nfyp,nfypyp,nfypy,nfypxp,nfypx,nfyyp,nfyy,nfyxp,nfyx,nfxpyp,nfxpy,nfxpxp,nfxpx,nfxyp,nfxy,nfxxp,nfxx,hx,gx,gxx,eta);

x0 = [0 0];
[Y, X] = simu_2nd(gx, hx, gxx, hxx, gss, hss, eta, sig, x0, e);

smpl = 0.9*num_sim:num_sim;

figure(10); plot(smpl,X(smpl,1),smpl,Y(smpl,3),smpl,Y(smpl,4),smpl,Y(smpl,5),smpl,Y(smpl,6));
legend('Capital','Consumption','Labor','Investment','Output')