neoclassical_model_run;

sig = 0.01;
num_sim = 1000;
e = randn(num_sim, 1);


x0 = [k a];
[Y, X] = simu_2nd(gx, hx, gxx, hxx, gss, hss, eta, sig, x0, e);
figure(1); plot(X(:,1)); title('Capital')
figure(5); plot(Y); title('Consumption')