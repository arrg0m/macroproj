[fx,fxp,fy,fyp,fypyp,fypy,fypxp,fypx,fyyp,fyy,fyxp,fyx,fxpyp,fxpy,fxpxp,fxpx,fxyp,fxy,fxxp,fxx,f] = EZGmodel(isLog);


%Numerical Evaluation
%Steady State and Parameter Values
[par_beta, par_delta, par_gamma, par_alpha, par_phi,...
    pf, z, w, bs, bn, n, k, invt, out, cs, cn, Js, Rs, Jn, Rn,...
    pfp, zp, wp, bsp, bnp, np, kp, invtp, outp, csp, cnp, Jsp, Rsp, Jnp, Rnp] = EZGmodel_ss(isLog, par_lambda);
%Order of approximation desired 
approx = 2;
eta = [0 1]';

%Obtain numerical derivatives of f
num_eval;

%First-order approximation
[gx,hx] = gx_hx(nfy,nfx,nfyp,nfxp);

%Second-order approximation
[gxx,hxx] = gxx_hxx(nfx,nfxp,nfy,nfyp,nfypyp,nfypy,nfypxp,nfypx,nfyyp,nfyy,nfyxp,nfyx,nfxpyp,nfxpy,nfxpxp,nfxpx,nfxyp,nfxy,nfxxp,nfxx,hx,gx);

[gss,hss] = gss_hss(nfx,nfxp,nfy,nfyp,nfypyp,nfypy,nfypxp,nfypx,nfyyp,nfyy,nfyxp,nfyx,nfxpyp,nfxpy,nfxpxp,nfxpx,nfxyp,nfxy,nfxxp,nfxx,hx,gx,gxx,eta);