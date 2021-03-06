isLog = 1;
par_lambda = 0.5; % 0.5, 1, 2
par_rho = 1/1.5;

EZGmodel_run;

sig = 0.01;
num_sim = 10000;
e = randn(num_sim, 1);

x0 = [0 0];
[Y, X] = simu_2nd(gx, hx, gxx, hxx, gss, hss, eta, sig, x0, e);

smpl = 0.9*num_sim:num_sim;

figure(10); plot(smpl,X(smpl,1),smpl,Y(smpl,5),smpl,Y(smpl,6),smpl,Y(smpl,9),'--',smpl,Y(smpl,10));
legend('Capital','Consumption (Shareholders)','Consumption (Workers)','Investment','Output')