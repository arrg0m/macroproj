
[fx,fxp,fy,fyp,fypyp,fypy,fypxp,fypx,fyyp,fyy,fyxp,fyx,fxpyp,fxpy,fxpxp,fxpx,fxyp,fxy,fxxp,fxx,f] = EZRBCmodel(isLog);

[par_beta, par_delta, par_gamma, par_alpha, par_phi, eta,...
            con, invt, out, k, z, n, J, R, ...
            conp, invtp, outp, kp, zp, np, Jp, Rp] = EZRBCmodel_ss(isLog, par_lambda);
approx = 2;
num_eval;

[gx,hx] = gx_hx(nfy,nfx,nfyp,nfxp);
[gxx,hxx] = gxx_hxx(nfx,nfxp,nfy,nfyp,nfypyp,nfypy,nfypxp,nfypx,nfyyp,nfyy,nfyxp,nfyx,nfxpyp,nfxpy,nfxpxp,nfxpx,nfxyp,nfxy,nfxxp,nfxx,hx,gx);
[gss,hss] = gss_hss(nfx,nfxp,nfy,nfyp,nfypyp,nfypy,nfypxp,nfypx,nfyyp,nfyy,nfyxp,nfyx,nfxpyp,nfxpy,nfxpxp,nfxpx,nfxyp,nfxy,nfxxp,nfxx,hx,gx,gxx,eta);