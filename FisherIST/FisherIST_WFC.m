% Two-shock case (with initialization)
FisherIST_run;
gss_av = gss(5);

% Neutral shock only
eta = [0 0; 1-par_alpha 0; 0 0];
[gss,hss] = gss_hss(nfx,nfxp,nfy,nfyp,nfypyp,nfypy,nfypxp,nfypx,nfyyp,nfyy,nfyxp,nfyx,nfxpyp,nfxpy,nfxpxp,nfxpx,nfxyp,nfxy,nfxxp,nfxx,hx,gx,gxx,eta);
gss_a = gss(5);

% INVT shock only
eta = [0 0; 0 0; 0 1];
[gss,hss] = gss_hss(nfx,nfxp,nfy,nfyp,nfypyp,nfypy,nfypxp,nfypx,nfyyp,nfyy,nfyxp,nfyx,nfxpyp,nfxpy,nfxpxp,nfxpx,nfxyp,nfxy,nfxxp,nfxx,hx,gx,gxx,eta);
gss_v = gss(5);

% Model-assumed variances
sig_v = 0.01;
sig_a = 0.01*(1-par_alpha);

% Welfare cost evaluation: direct
lambda_a  = 1-exp(0.5*(1-par_beta)*gss_a*sig_a^2);
lambda_v  = 1-exp(0.5*(1-par_beta)*gss_v*sig_v^2);
lambda_av = 1-exp(0.5*(1-par_beta)*gss_av*(sig_a^2+sig_v^2));

% Welfare cost: indirect (`conditional')
lambda_v_a = 1-(1-lambda_av)/(1-lambda_a);
lambda_a_v = 1-(1-lambda_av)/(1-lambda_v);

% Welfare costs in percentage units
lambda_a_perc  = 100*lambda_a;
lambda_v_perc  = 100*lambda_v;
lambda_av_perc = 100*lambda_av;
lambda_v_a_perc = 100*lambda_v_a;
lambda_a_v_perc = 100*lambda_a_v;