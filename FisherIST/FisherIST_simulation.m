FisherIST_run;

num_sim = 100;
Uav = zeros(num_sim, 1);
Ua = zeros(num_sim, 1);
Uv = zeros(num_sim, 1);
U = zeros(num_sim, 1);

for itr_sim = 1:num_sim;
    sig = 0.01;
    num_sim = 1000;
    e = randn(num_sim, 2); % When both shock exists

    x0 = [0, 0, 0];
    [Yav, Xav] = simu_2nd(gx, hx, gxx, hxx, gss, hss, eta, sig, x0, e);
    lnA = par_gamma*(1:num_sim)' + sig*cumsum(e(:,1));
    lnV = par_nu*(1:num_sim)' + sig*cumsum(e(:,2));
    lnQ = (1/(1-par_alpha))*lnA + (par_alpha/(1-par_alpha))*lnV;
    % lnX = Yav(2:end, 3) + lnQ; 
    lnC = Yav(2:end, 1) + lnQ; 
    H = Yav(2:end, 2); 
    % lnY = Yav(2:end, 4) + lnQ;
    Uav(itr_sim) = sum((lnC - H) .* par_beta.^((1:num_sim)-1)');

    ee = e;
    ee(:, 2) = 0;
    [Ya, Xa] = simu_2nd(gx, hx, gxx, hxx, gss, hss, eta, sig, x0, ee);
    lnA = par_gamma*(1:num_sim)' + sig*cumsum(e(:,1));
    lnV = par_nu*(1:num_sim)' + sig*cumsum(e(:,2));
    lnQ = (1/(1-par_alpha))*lnA + (par_alpha/(1-par_alpha))*lnV;
    lnC = Ya(2:end, 1) + lnQ; 
    H = Ya(2:end, 2); 
    Ua(itr_sim) = sum((lnC - H) .* par_beta.^((1:num_sim)-1)');

    ee = e;
    ee(:, 1) = 0;
    [Yv, Xv] = simu_2nd(gx, hx, gxx, hxx, gss, hss, eta, sig, x0, ee);
    lnA = par_gamma*(1:num_sim)' + sig*cumsum(e(:,1));
    lnV = par_nu*(1:num_sim)' + sig*cumsum(e(:,2));
    lnQ = (1/(1-par_alpha))*lnA + (par_alpha/(1-par_alpha))*lnV;
    lnC = Yv(2:end, 1) + lnQ; 
    H = Yv(2:end, 2); 
    Uv(itr_sim) = sum((lnC - H) .* par_beta.^((1:num_sim)-1)');

    ee = zeros(num_sim, 2);
    [Y, X] = simu_2nd(gx, hx, gxx, hxx, gss, hss, eta, sig, x0, ee);
    lnA = par_gamma*(1:num_sim)' + sig*cumsum(e(:,1));
    lnV = par_nu*(1:num_sim)' + sig*cumsum(e(:,2));
    lnQ = (1/(1-par_alpha))*lnA + (par_alpha/(1-par_alpha))*lnV;
    lnC = Y(2:end, 1) + lnQ; 
    H = Y(2:end, 2); 
    U(itr_sim) = sum((lnC - H) .* par_beta.^((1:num_sim)-1)');
end

uav = mean(Uav);
ua = mean(Ua);
uv = mean(Uv);
u = mean(U);

g_hat = par_alpha*par_nu/(1-par_beta);
d_hat = 0.5*(par_alpha/(1-par_alpha))^2 * 0.01^2;

lambda_invt = exp((1-par_beta)*(d_hat/(1-par_beta) + ua - uav)) - 1;