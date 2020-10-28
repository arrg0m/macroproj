% generalize codes?
% 1. number of variables, to detrend or not for each, ...
% 2. number of lags
% 3. date structures
% 4. IRF?
% ...nope

% suptitle function?


% flags
useunitvec = 1; % for non-differenced variables...?
source = 0;

%% BLOCK00 IMPORT
if source == 1
    lambdaIV = 10^(-9);
    % US Data: Date from 1947Q2 to 2015Q4
    FisherIV_import;
    year_begin = 1947;
    quat_begin = 2;
    %begin_smpl = [1955, 1]; end_smpl = [1979, 2];
    
    %begin_smpl = [1983, 2]; end_smpl = [2000, 4];
    begin_smpl = [2001, 1]; end_smpl = [2014, 4];
    diff_h = 0;
    diff_r = 0;
    diff_pi = 0;
elseif source == 0
%    lambdaIV = 10^(-12);
    [tmp1, tmp2, raw] = xlsread('C:\Users\hmbyun\Dropbox\Work\ishock_Data_KR.xlsx','Export','J2:O63');
    R = cellfun(@(x) ~isnumeric(x) || isnan(x),raw); % Find non-numeric cells
    raw(R) = {0.0}; % Replace non-numeric cells
    data = cell2mat(raw);
    dp = data(:,1);
    da = data(:,2);
    h = data(:,3);
    r = data(:,4);
    pi = data(:,5);
    re = data(:,6); % re-specification?
	clearvars data raw R columnIndices;
    
    year_begin = 2002;
    quat_begin = 1;
    begin_smpl = [2002, 2]; end_smpl = [2016, 4];
    % seasonality issue?
    % 일단 bloomberg data 써보자 (if available...)

%     %Option 02 (non-differencing results in large seasonality)
%     lambdaIV = 10^(-7);
%     diff_h = 0;
%     diff_r = 0;
%     diff_pi = 0;

      lambdaIV = 10^(-16);
      diff_h = 1;
      diff_r = 0;
      diff_pi = 1;

      r = re;
      %Option 01 (Presentation)
%       lambdaIV = (1/10)^4;
%       diff_h = 0;
%       diff_r = 1;
%       diff_pi = 1;

end

% dummy
len = size(dp,1);
datevec = zeros(len, 2);
datevec(1,:) = [year_begin, quat_begin];
year = year_begin;
quat = quat_begin;
begin_itr = 1;
end_itr = 1;
for itr = 2:len
    if quat == 4
        quat = 1;
        year = year + 1;
    else
        quat = quat + 1;
    end
    datevec(itr,:) = [year, quat];
    
    % possible to be optimized (via direct search)
    if year == begin_smpl(1) && quat == begin_smpl(2)
        begin_itr = itr;
    end
    if year == end_smpl(1) && quat == end_smpl(2)
        end_itr = itr;
    end
end

smpl = begin_itr:end_itr;

%differencing h
%h = h - [0; h(1:end-1)];
%이거 차분하면 irf도 바꿔야함

% differencing r; or detrend? HP-filtered residual??
if diff_h == 1
    h = h - [0; h(1:end-1)];
end
if diff_r == 1
    r = r - [0; r(1:end-1)];
end

% differencing pi; or detrend?
if diff_pi == 1
    pi = pi - [0; pi(1:end-1)];
end

is_diff5 = [1 1 diff_h diff_r];
is_diff4 = [1 diff_h diff_r];
% Q) Is demean necessary? - add columns with ones instead (hope it works)
% dp = dp - mean(dp(smpl));
% da = da - mean(da(smpl));
% h  = (h  - mean(h(smpl)))/std(h(smpl));
% h = h - mean(h(smpl));
% r  = r  - mean(r(smpl) );
% pi = pi - mean(pi(smpl));
% figure(1);
% subplot(221);
% plot(smpl, dp(smpl));
% ylabel('dp');
% subplot(222);
% plot(smpl, da(smpl));
% ylabel('da');
% subplot(223);
% plot(smpl, h(smpl));
% ylabel('h')
% subplot(224);
% plot(smpl, r(smpl));
% ylabel('r')

% data range problem...?

% ...r, pi exogenous?

%% BLOCK01 PREPROCESSING (DERIVED VARIABLES)
dp1 = [0; dp(1:end-1)];
dp2 = [0; 0; dp(1:end-2)]; % or dp2 = [0; dp1(1:end-1)];
dp3 = [0; 0; 0; dp(1:end-3)];
dp4 = [0; 0; 0; 0; dp(1:end-4)];
ddp = dp - dp1;
ddp1 = dp1 - dp2;
ddp2 = dp2 - dp3;
ddp3 = dp3 - dp4;

da1 = [0; da(1:end-1)];
da2 = [0; 0; da(1:end-2)];
da3 = [0; 0; 0; da(1:end-3)];
da4 = [0; 0; 0; 0; da(1:end-4)];
dda = da - da1;
dda1 = da1 - da2;
dda2 = da2 - da3;
dda3 = da3 - da4;

h1 = [0; h(1:end-1)];
h2 = [0; 0; h(1:end-2)];
h3 = [0; 0; 0; h(1:end-3)];
h4 = [0; 0; 0; 0; h(1:end-4)];
dh = h - h1;
dh1 = h1 - h2;
dh2 = h2 - h3;
dh3 = h3 - h4;

r1 = [0; r(1:end-1)];
r2 = [0; 0; r(1:end-2)];
r3 = [0; 0; 0; r(1:end-3)];
r4 = [0; 0; 0; 0; r(1:end-4)];
dr = r - r1;
dr1 = r1 - r2;
dr2 = r2 - r3;
dr3 = r3 - r4;

if useunitvec == 1
    unitvec = ones(len, 1);
end
%% BLOCK02 INSTRUMENTAL VARIABLE ESTIMATION
Z0 = [dp1, dp2, dp3, dp4, da1, da2, da3, da4, h1, h2, h3, h4, r1, r2, r3, r4];

% Equation 1: on dp
y = dp; y = y(smpl);
X = [dp1, dp2, dp3, dp4, dda, dda1, dda2, dda3, dh, dh1, dh2, dh3, dr, dr1, dr2, dr3];
if useunitvec == 1
    X = [X, unitvec];
end
X = X(smpl,:);
Z = Z0(smpl,:);
[coef1, resid1] = myIV(y, X, Z, lambdaIV);

% Equation 2: on da
y = da; y = y(smpl);
X = [dp, dp1, dp2, dp3, dp4, da1, da2, da3, da4, dh, dh1, dh2, dh3, dr, dr1, dr2, dr3]; 
if useunitvec == 1
    X = [X, unitvec];
end
X = X(smpl,:);
Z = [Z0(smpl, :), resid1];
[coef2, resid2] = myIV(y, X, Z, lambdaIV);

% Equation 3: on h
y = h; y = y(smpl);
X = [dp, dp1, dp2, dp3, da, da1, da2, da3, h1, h2, h3, r, r1, r2, r3]; 
if useunitvec == 1
    X = [X, unitvec];
end
X = X(smpl,:);
Z = [Z0(smpl,:), resid1, resid2]; 
[coef3, resid3] = myIV(y, X, Z, lambdaIV);

% Equation 4: on r
y = r; y = y(smpl);
X = [dp, dp1, dp2, dp3, da, da1, da2, da3, h, h1, h2, h3, r1, r2, r3];
if useunitvec == 1
    X = [X, unitvec];
end
X = X(smpl,:);
Z = [Z0(smpl,:), resid1, resid2]; 
[coef4, resid4] = myIV(y, X, Z, lambdaIV); 

% Coefficient Recovery
A5 = [1, -coef1(5), -coef1(9), -coef1(13);
     -coef2(1), 1, -coef2(10), -coef2(14);
     -coef3(1), -coef3(5), 1, -coef3(12);
     -coef4(1), -coef4(5), -coef4(9), 1];

G51 = [coef1(1), -coef1(5) + coef1(6), -coef1(9) + coef1(10), -coef1(13) + coef1(14);
        coef2(2), coef2(6), -coef2(10) + coef2(11), -coef2(14) + coef2(15);
        -coef3(1) + coef3(2), -coef3(5) + coef3(6),  1 + coef3(9), -coef3(12) + coef3(13);
        -coef4(1) + coef4(2), -coef4(5) + coef4(6), -coef4(9) + coef4(10), 1 + coef4(13)];
    
G52 = [coef1(2), -coef1(6) + coef1(7), -coef1(10) + coef1(11), -coef1(14) + coef1(15);
        coef2(3), coef2(7), -coef2(11) + coef2(12), -coef2(15) + coef2(16);
        -coef3(2) + coef3(3), -coef3(6) + coef3(7), - coef3(9) + coef3(10), -coef3(13) + coef3(14);
        -coef4(2) + coef4(3), -coef4(6) + coef4(7), -coef4(10) + coef4(11), -coef4(13) + coef4(14)];

G53 = [coef1(3), -coef1(7) + coef1(8), -coef1(11) + coef1(12), -coef1(15) + coef1(16);
        coef2(4), coef2(8), -coef2(12) + coef2(13), -coef2(16) + coef2(17);
        -coef3(3) + coef3(4), -coef3(7) + coef3(8), -coef3(10) + coef3(11), -coef3(14) + coef3(15);
        -coef4(3) + coef4(4), -coef4(7) + coef4(8), -coef4(11) + coef4(12), -coef4(14) + coef4(15)];

G54 = [coef1(4), -coef1(8), -coef1(12), -coef1(16);
        coef2(5), coef2(9), -coef2(13), -coef2(17);
        -coef3(4), -coef3(8), -coef3(11), -coef3(15);
        -coef4(4), -coef4(8), -coef4(12), -coef4(15)];

B5sq = cov(resid1, resid2); B5 = B5sq^0.5;
Sig5 = cov([resid1, resid2, resid3, resid4]);

cn5 = zeros(4, 1);
if useunitvec == 1
    cn5 = [coef1(17); coef2(18); coef3(16); coef4(16)];
end

%% BLOCK03 THEORETICAL IRF COMPUTATION
C5 = zeros(4, 4, 48); % coordinates: to, from, period
C5(:, :, 1) = inv(A5);
C5(:, :, 2) = C5(:, :, 1)*(G51*C5(:, :, 1));
C5(:, :, 3) = C5(:, :, 1)*(G52*C5(:, :, 1) + G51*C5(:, :, 2));
C5(:, :, 4) = C5(:, :, 1)*(G53*C5(:, :, 1) + G52*C5(:, :, 2) + G51*C5(:, :, 3));
C5(:, :, 5) = C5(:, :, 1)*(G54*C5(:, :, 1) + G53*C5(:, :, 2) + G52*C5(:, :, 3) + G51*C5(:, :, 4));
for itr = 6:48
    C5(:, :, itr) = C5(:, :, 1)*(G54*C5(:, :, itr-4) + G53*C5(:, :, itr-3) + G52*C5(:, :, itr-2) + G51*C5(:, :, itr-1));
end

for itr = 1:48
    C5(:, :, itr) = C5(:, :, itr)*Sig5^0.5;
end
FEVD5 = myFEVD(C5, is_diff5);


%% BLOCK02' INSTRUMENTAL VARIABLE ESTIMATION WITHOUT DP
Z0 = [da1, da2, da3, da4, h1, h2, h3, h4, r1, r2, r3, r4];

% Equation 2: on da
y = da; y = y(smpl);
X = [da1, da2, da3, da4, dh, dh1, dh2, dh3, dr, dr1, dr2, dr3];
if useunitvec == 1
    X = [X, unitvec];
end
X = X(smpl,:);
Z = Z0(smpl, :);
[coef2, resid2] = myIV(y, X, Z, lambdaIV);

% Equation 3: on h
y = h; y = y(smpl);
X = [da, da1, da2, da3, h1, h2, h3, r, r1, r2, r3];
if useunitvec == 1
    X = [X, unitvec];
end
X = X(smpl,:);
Z = [Z0(smpl,:), resid2]; 
[coef3, resid3] = myIV(y, X, Z, lambdaIV);

% Equation 4: on r
y = r; y = y(smpl);
X = [da, da1, da2, da3, h, h1, h2, h3, r1, r2, r3,];
if useunitvec == 1
    X = [X, unitvec];
end
X = X(smpl,:);
Z = [Z0(smpl,:), resid2]; 
[coef4, resid4] = myIV(y, X, Z, lambdaIV); 


A4 = [1, -coef2(5), -coef2(9); 
    -coef3(1), 1, -coef3(8);
    -coef4(1), -coef4(5), 1];

G41 = [coef2(1), -coef2(5) + coef2(6), -coef2(9) + coef2(10); 
    -coef3(1) + coef3(2), 1 + coef3(5), -coef3(8) + coef3(9);
    -coef4(1) + coef4(2), -coef4(5) + coef4(6), 1 + coef4(9)];

G42 = [coef2(2), -coef2(6) + coef2(7), -coef2(10) + coef2(11); 
    -coef3(2) + coef3(3), -coef3(5) + coef3(6), -coef3(9) + coef3(10);
    -coef4(2) + coef4(3), -coef4(6) + coef4(7), -coef4(9) + coef4(10)];

G43 = [coef2(3), -coef2(7) + coef2(8), -coef2(11) + coef2(12); 
    -coef3(3) + coef3(4), -coef3(6) + coef3(7), -coef3(10) + coef3(11);
    -coef4(3) + coef4(4), -coef4(7) + coef4(8), -coef4(10) + coef4(11)];

G44 = [coef2(4), -coef2(8), -coef2(12); 
    -coef3(4), -coef3(7), -coef3(11);
    -coef4(4), -coef4(8), -coef4(11)];

B4sq = var(resid2); B4 = B4sq^0.5;
Sig4 = cov([resid2, resid3, resid4]);

cn4 = [0; 0; 0];
if useunitvec == 1
    cn4 = [coef2(13); coef3(12); coef4(12)];
end
%% BLOCK03' THEORETICAL IRF COMPUTATION
C4 = zeros(3, 3, 48);
C4(:, :, 1) = inv(A4);
C4(:, :, 2) = A4\G41/A4;
C4(:, :, 3) = A4\(G42/A4 + G41*C4(:, :, 2));
C4(:, :, 4) = A4\(G43/A4 + G42*C4(:, :, 2) + G41*C4(:, :, 3));
C4(:, :, 5) = A4\(G44/A4 + G43*C4(:, :, 2) + G42*C4(:, :, 3) + G41*C4(:, :, 4));
for itr = 6:48
    C4(:, :, itr) = A4\(G44*C4(:, :, itr-4) + G43*C4(:, :, itr-3) + G42*C4(:, :, itr-2) + G41*C4(:, :, itr-1));
end

for itr = 1:48
    C4(:, :, itr) = C4(:, :, itr)*Sig4^0.5;
end
FEVD4 = myFEVD(C4, is_diff4);

% %% GRAPH
% figure(2);
% 
% subplot(221)
% dp_dp = squeeze(C5(1,1,:))*B5(1,1);
% p_dp = cumsum(dp_dp);
% dp_da = squeeze(C5(1,2,:))*B5(2,2);
% p_da = cumsum(dp_da);
% plot(0:48, [0; -p_dp], 'b-', 0:48, [0; p_da], 'r--')
% ylabel('Response of investment price')
% legend('V', 'A')
% 
% da_dp = squeeze(C5(2,1,:))*B5(1,1);
% a_dp = cumsum(da_dp);
% da_da = squeeze(C5(2,2,:))*B5(2,2);
% a_da = cumsum(da_da);
% da_daonly = squeeze(C4(1,1,:))*B4;
% a_daonly = cumsum(da_daonly);
% subplot(222)
% plot(0:48, [0; -a_dp], 'b-', 0:48, [0; a_da], 'r--', 0:48, [0; a_daonly], 'k:')
% ylabel('Response of productivity')
% legend('V', 'A', 'A-only')
% 
% h_dp = squeeze(C5(3,1,:))*B5(1,1);
% h_da = squeeze(C5(3,2,:))*B5(2,2);
% h_daonly = squeeze(C4(2,1,:))*B4;
% if diff_h == 1
%     dh_dp = h_dp;
%     dh_da = h_da;
%     dh_daonly = h_daonly;
%     h_dp = cumsum(dh_dp);
%     h_da = cumsum(dh_da);
%     h_daonly = cumsum(dh_daonly);
% end
% subplot(223)
% plot(0:48, [0; -h_dp], 'b-', 0:48, [0; h_da], 'r--', 0:48, [0; h_daonly], 'k:')
% ylabel('Response of Hours worked')
% legend('V', 'A', 'A-only')
% 
% y_dp = a_dp + h_dp;
% y_da = a_da + h_da;
% y_daonly = a_daonly + h_daonly;
% subplot(224)
% plot(0:48, [0; -y_dp], 'b-', 0:48, [0; y_da], 'r--', 0:48, [0; y_daonly], 'k:')
% ylabel('Response of Output')
% legend('V', 'A', 'A-only')
%suptitle('VMA Coefficient-based IRF')

% TODO FEVD

%% BLOCK03 VER. BA

% B = diag([std(resid1), std(resid2), std(resid3), std(resid4), std(resid5)]);
% %Bsq = cov([resid1, resid2, resid3, resid4, resid5]); B = Bsq^0.5;
% 
% BA = B\A;
% BG1 = B\G1;
% BG2 = B\G2;
% BG3 = B\G3;
% BG4 = B\G4;

% C0 = inv(BA);
% C = zeros(5, 5, 48);
% C(:, :, 1) = C0*BG1*C0;
% C(:, :, 2) = C0*(BG2*C0 + BG1*C(:, :, 1));
% C(:, :, 3) = C0*(BG3*C0 + BG2*C(:, :, 1) + BG1*C(:, :, 2));
% C(:, :, 4) = C0*(BG4*C0 + BG3*C(:, :, 1) + BG2*C(:, :, 2) + BG1*C(:, :, 3));
% for itr = 5:48
%     C(:, :, itr) = C0*(BG4*C(:, :, itr-4) + BG3*C(:, :, itr-3) + BG2*C(:, :, itr-2) + BG1*C(:, :, itr-1));
% end
% 
% dp_dp = squeeze(C(1,1,:));
% p_dp = cumsum(dp_dp);
% dp_da = squeeze(C(1,2,:));
% p_da = cumsum(dp_da);
% plot(1:48, -p_dp, 1:48, p_da)
% da_dp = squeeze(C(2,1,:));
% a_dp = cumsum(da_dp);
% da_da = squeeze(C(2,2,:));
% a_da = cumsum(da_da);
% plot(1:48, -a_dp, 1:48, a_da)
% h_dp = squeeze(C(3,1,:));
% h_da = squeeze(C(3,2,:));
% plot(1:48, -h_dp, 1:48, h_da)

%% BLOCK03' SIMULATED IRF COMPUTATION
% http://www.dynare.org/DynareWiki/IrFs

num_sim = 100;
len_pre = 100;
len_irf = 50;
len_sim = len_pre + len_irf;
smpl_irf = (len_pre+1):len_sim;
smpl_irf_plt = 1:len_irf;

% 5-var

y0 = zeros(4, 1);
y = zeros(4, len_sim, num_sim);
y1 = zeros(4, len_sim, num_sim);
y2 = zeros(4, len_sim, num_sim);

for itr = 1:num_sim
    e = Sig5^0.5 * randn(4, len_sim);
    e1 = e;
    e2 = e;
    e1(1, len_pre+1) = e1(1, len_pre+1) - Sig5(1,1);
    e2(2, len_pre+1) = e2(2, len_pre+1) + Sig5(2,2);
    y(:, 1, itr) = A5\cn5 + A5\(G51*y0) + A5\e(:, 1);
    y(:, 2, itr) = A5\cn5 + A5\(G51*y(:, 1, itr) + G52*y0) + A5\e(:, 2);
    y(:, 3, itr) = A5\cn5 + A5\(G51*y(:, 2, itr) + G52*y(:, 1, itr) + G53*y0) + A5\e(:,3);
    y(:, 4, itr) = A5\cn5 + A5\(G51*y(:, 3, itr) + G52*y(:, 2, itr) + G53*y(:, 1, itr) + G54*y0) + A5\e(:,4);
    for itr2 = 5:len_sim
        y(:, itr2, itr) = A5\cn5 + A5\(G51*y(:, itr2-1, itr) + G52*y(:, itr2-2, itr) + G53*y(:, itr2-3, itr) + G54*y(:, itr2-4, itr)) + A5\e(:,itr2);
    end
    
    y1(:, 1, itr) = A5\cn5 + A5\(G51*y0) + A5\e1(:, 1);
    y1(:, 2, itr) = A5\cn5 + A5\(G51*y1(:, 1, itr) + G52*y0) + A5\e1(:, 2);
    y1(:, 3, itr) = A5\cn5 + A5\(G51*y1(:, 2, itr) + G52*y1(:, 1, itr) + G53*y0) + A5\e1(:,3);
    y1(:, 4, itr) = A5\cn5 + A5\(G51*y1(:, 3, itr) + G52*y1(:, 2, itr) + G53*y1(:, 1, itr) + G54*y0) + A5\e1(:,4);
    for itr2 = 5:len_sim
        y1(:, itr2, itr) = A5\cn5 + A5\(G51*y1(:, itr2-1, itr) + G52*y1(:, itr2-2, itr) + G53*y1(:, itr2-3, itr) + G54*y1(:, itr2-4, itr)) + A5\e1(:,itr2);
    end   
    
    y2(:, 1, itr) = A5\cn5 + A5\(G51*y0) + A5\e2(:, 1);
    y2(:, 2, itr) = A5\cn5 + A5\(G51*y2(:, 1, itr) + G52*y0) + A5\e2(:, 2);
    y2(:, 3, itr) = A5\cn5 + A5\(G51*y2(:, 2, itr) + G52*y2(:, 1, itr) + G53*y0) + A5\e2(:,3);
    y2(:, 4, itr) = A5\cn5 + A5\(G51*y2(:, 3, itr) + G52*y2(:, 2, itr) + G53*y2(:, 1, itr) + G54*y0) + A5\e2(:,4);
    for itr2 = 5:len_sim
        y2(:, itr2, itr) = A5\cn5 + A5\(G51*y2(:, itr2-1, itr) + G52*y2(:, itr2-2, itr) + G53*y2(:, itr2-3, itr) + G54*y2(:, itr2-4, itr)) + A5\e2(:,itr2);
    end
end
y_res = mean(y, 3);
y1_res = mean(y1, 3);
y2_res = mean(y2, 3);



% 4-var
z0 = zeros(3, 1);
z = zeros(3, len_sim, num_sim);
z1 = zeros(3, len_sim, num_sim);

for itr = 1:num_sim
    f = Sig4^0.5 * randn(3, len_sim);
    f1 = f;
    f1(1, len_pre+1) = f1(1, len_pre+1) + Sig4(1,1);
    z(:, 1, itr) = A4\cn4 + A4\(G41*z0) + A4\f(:, 1);
    z(:, 2, itr) = A4\cn4 + A4\(G41*z(:, 1, itr) + G42*z0) + A4\f(:, 2);
    z(:, 3, itr) = A4\cn4 + A4\(G41*z(:, 2, itr) + G42*z(:, 1, itr) + G43*z0) + A4\f(:,3);
    z(:, 4, itr) = A4\cn4 + A4\(G41*z(:, 3, itr) + G42*z(:, 2, itr) + G43*z(:, 1, itr) + G44*z0) + A4\f(:,4);
    for itr2 = 5:len_sim
        z(:, itr2, itr) = A4\cn4 + A4\(G41*z(:, itr2-1, itr) + G42*z(:, itr2-2, itr) + G43*z(:, itr2-3, itr) + G44*z(:, itr2-4, itr)) + A4\f(:,itr2);
    end
    
    z1(:, 1, itr) = A4\cn4 + A4\(G41*z0) + A4\f1(:, 1);
    z1(:, 2, itr) = A4\cn4 + A4\(G41*z1(:, 1, itr) + G42*z0) + A4\f1(:, 2);
    z1(:, 3, itr) = A4\cn4 + A4\(G41*z1(:, 2, itr) + G42*z1(:, 1, itr) + G43*z0) + A4\f1(:,3);
    z1(:, 4, itr) = A4\cn4 + A4\(G41*z1(:, 3, itr) + G42*z1(:, 2, itr) + G43*z1(:, 1, itr) + G44*z0) + A4\f1(:,4);
    for itr2 = 5:len_sim
        z1(:, itr2, itr) = A4\cn4 + A4\(G41*z1(:, itr2-1, itr) + G42*z1(:, itr2-2, itr) + G43*z1(:, itr2-3, itr) + G44*z1(:, itr2-4, itr)) + A4\f1(:,itr2);
    end
end
z_res = mean(z, 3);
z1_res = mean(z1, 3);

% collect

myirf = zeros(3, 4, len_irf);
myirf(1, 1, :) = cumsum(y1_res(1,smpl_irf) - y_res(1,smpl_irf)); % investment price response
myirf(2, 1, :) = cumsum(y2_res(1,smpl_irf) - y_res(1,smpl_irf));
myirf(3, 1, :) = 0;
myirf(1, 2, :) = cumsum(y1_res(2,smpl_irf) - y_res(2,smpl_irf)); % productivity response
myirf(2, 2, :) = cumsum(y2_res(2,smpl_irf) - y_res(2,smpl_irf));
myirf(3, 2, :) = cumsum(z1_res(1,smpl_irf) - z_res(1,smpl_irf));
myirf(1, 3, :) = y1_res(3,smpl_irf) - y_res(3,smpl_irf); % labor hour response
myirf(2, 3, :) = y2_res(3,smpl_irf) - y_res(3,smpl_irf);
myirf(3, 3, :) = z1_res(2,smpl_irf) - z_res(2,smpl_irf);
if diff_h == 1
    myirf(1, 3, :) = cumsum(myirf(1, 3, :));
    myirf(2, 3, :) = cumsum(myirf(2, 3, :));
    myirf(3, 3, :) = cumsum(myirf(3, 3, :));
end
myirf(1, 4, :) = myirf(1, 2, :) + myirf(1, 3, :); % output = productivity * hours
myirf(2, 4, :) = myirf(2, 2, :) + myirf(2, 3, :);
myirf(3, 4, :) = myirf(3, 2, :) + myirf(3, 3, :);
myirf = myirf*100; % convert to percentage unit

figure(3);
subplot(221);
plot(smpl_irf_plt, squeeze(myirf(1,1,:)), 'b-', smpl_irf_plt, squeeze(myirf(2,1,:)), 'r--', smpl_irf_plt, squeeze(myirf(3,1,:)), 'k:');
ylabel('Response of investment price')
legend('V', 'A', 'A-only')
subplot(222);
plot(smpl_irf_plt, squeeze(myirf(1,2,:)), 'b-', smpl_irf_plt, squeeze(myirf(2,2,:)), 'r--', smpl_irf_plt, squeeze(myirf(3,2,:)), 'k:');
ylabel('Response of productivity') % 이거 굳이 그려야 하나?
legend('V', 'A', 'A-only')
subplot(223);
plot(smpl_irf_plt, squeeze(myirf(1,3,:)), 'b-', smpl_irf_plt, squeeze(myirf(2,3,:)), 'r--', smpl_irf_plt, squeeze(myirf(3,3,:)), 'k:');
ylabel('Response of Hours worked')
legend('V', 'A', 'A-only')
subplot(224);
plot(smpl_irf_plt, squeeze(myirf(1,4,:)), 'b-', smpl_irf_plt, squeeze(myirf(2,4,:)), 'r--', smpl_irf_plt, squeeze(myirf(3,4,:)), 'k:');
ylabel('Response of Output')
legend('V', 'A', 'A-only')
%suptitle('Simulation-based IRF')


figure(4);
plot(1:48, squeeze(FEVD5(2, 1, :))', 'b-', 1:48, squeeze(FEVD5(2, 2, :))', 'r--', 1:48, squeeze(FEVD4(1, 1, :))', 'k:')
title('Forecast Error Variance Decomposition: Productivity')
legend('V', 'A', 'A-only')

figure(5);
plot(1:48, squeeze(FEVD5(3, 1, :))', 'b-', 1:48, squeeze(FEVD5(3, 2, :))', 'r--', 1:48, squeeze(FEVD4(2, 1, :))', 'k:')
title('Forecast Error Variance Decomposition: Hours Worked')
legend('V', 'A', 'A-only')


% Well, I think this is quite more accurate, due to the presence of the trends
% (constant terms)


% %% BLOCK04 FEVD
% % New Intro., Section 2.3.3
% 
% % Five-variable case
% C5B = zeros(4, 4, 48); % j, k, t
% for itr = 1:48
%     C5B(:, :, itr) = C5(:, :, itr) * Sig5^0.5; % THETA from the textbook
% end
% C5B_contrib = zeros(4, 4, 48); % 2.3.36
% for itrj = 1:4
%     for itrk = 1:4
%         C5B_contrib(itrj, itrk, 1) = C5B(itrj, itrk, 1)^2;
%         for itrh = 2:48
%             C5B_contrib(itrj, itrk, itrh) = C5B_contrib(itrj, itrk, itrh-1) + C5B(itrj, itrk, itrh)^2;
%         end
%     end
% end
% MSE5B = zeros(4, 48);
% for itrj = 1:4
%     MSE5B(itrj, 1) = sum(C5B_contrib(itrj, :, 1));
%     for itrh = 2:48
%         MSE5B(itrj, itrh) = MSE5B(itrj, itrh-1) + sum(C5B_contrib(itrj, :, itrh));
%     end
% end
% C5B_fevd = zeros(4, 4, 48);
% for itrj = 1:4
%     for itrk = 1:4
%         for itrh = 1:48
%             C5B_fevd(itrj, itrk, itrh) = sum(C5B_contrib(itrj, itrk, 1:itrh))/MSE5B(itrj, itrh);
%         end
%     end
% end
% 
% % Four-variable Case
% C4B = zeros(3, 3, 48); % j, k, t
% for itr = 1:48
%     C4B(:, :, itr) = C4(:, :, itr)*Sig4^0.5;
% end
% C4B_contrib = zeros(3, 3, 48); % 2.3.36
% for itrj = 1:3
%     for itrk = 1:3
%         C4B_contrib(itrj, itrk, 1) = C4B(itrj, itrk, 1)^2;
%         for itrh = 2:48
%             C4B_contrib(itrj, itrk, itrh) = C4B_contrib(itrj, itrk, itrh-1) + C4B(itrj, itrk, itrh)^2;
%         end
%     end
% end
% MSE4B = zeros(3, 48);
% for itrj = 1:3
%     MSE4B(itrj, 1) = sum(C4B_contrib(itrj, :, 1));
%     for itrh = 2:48
%         MSE4B(itrj, itrh) = MSE4B(itrj, itrh-1) + sum(C4B_contrib(itrj, :, itrh));
%     end
% end
% C4B_fevd = zeros(3, 3, 48);
% for itrj = 1:3
%     for itrk = 1:3
%         for itrh = 1:48
%             C4B_fevd(itrj, itrk, itrh) = sum(C4B_contrib(itrj, itrk, 1:itrh))/MSE4B(itrj, itrh);
%         end
%     end
% end
% % 
% % figure(4);
% % plot(1:48, squeeze(C5B_fevd(2, 1, :))', 'b-', 1:48, squeeze(C5B_fevd(2, 2, :))', 'r--', 1:48, squeeze(C4B_fevd(1, 1, :))', 'k:')
% % title('Forecast Error Variance Decomposition: Productivity')
% % legend('V', 'A', 'A-only')
% % 
% % figure(5);
% % plot(1:48, squeeze(C5B_fevd(3, 1, :))', 'b-', 1:48, squeeze(C5B_fevd(3, 2, :))', 'r--', 1:48, squeeze(C4B_fevd(2, 1, :))', 'k:')
% % title('Forecast Error Variance Decomposition: Hours Worked')
% % legend('V', 'A', 'A-only')
% 
% 
% %%
% 
% MSE5_tot = zeros(4, 4, 48);
% MSE4_tot = zeros(3, 3, 48);
% MSE5_tot(:, :, 1) = C5B(:, :, 1) * C5B(:, :, 1)';
% MSE4_tot(:, :, 1) = C4B(:, :, 1) * C4B(:, :, 1)';
% for itrh = 2:48
%     MSE5_tot(:, :, itrh) = MSE5_tot(:, :, itrh-1) + C5B(:, :, itrh)*C5B(:, :, itrh)';
%     MSE4_tot(:, :, itrh) = MSE4_tot(:, :, itrh-1) + C4B(:, :, itrh)*C4B(:, :, itrh)';
% end
% 
% %%
% 
% X5 = zeros(4, 4, 48);
% X4 = zeros(3, 3, 48);
% for itrj = 1:4
%     for itrk = 1:4
%         X5(itrj, itrk, 1) = C5B(itrj, itrk, 1)^2;
%         for itrh = 2:48
%             X5(itrj, itrk, itrh) = X5(itrj, itrk, itrh-1) + C5B(itrj, itrk, itrh)^2;
%         end
%     end
% end
% for itrj = 1:3
%     for itrk = 1:3
%         X4(itrj, itrk, 1) = C4B(itrj, itrk, 1)^2;
%         for itrh = 2:48
%             X4(itrj, itrk, itrh) = X4(itrj, itrk, itrh-1) + C4B(itrj, itrk, itrh)^2;
%         end
%     end
% end
% 
% W5 = zeros(4, 4, 48);
% W4 = zeros(3, 3, 48);
% for itrj = 1:4
%     for itrh = 1:48
%         W5(itrj, :, itrh) = X5(itrj, :, itrh)/sum(X5(itrj, :, itrh));
%     end
% end
% for itrj = 1:3
%     for itrh = 1:48
%         W4(itrj, :, itrh) = X4(itrj, :, itrh)/sum(X4(itrj, :, itrh));
%     end
% end
% 
% 
% figure(4);
% plot(1:48, squeeze(W5(2, 1, :))', 'b-', 1:48, squeeze(W5(2, 2, :))', 'r--', 1:48, squeeze(W4(1, 1, :))', 'k:')
% title('Forecast Error Variance Decomposition: Productivity')
% legend('V', 'A', 'A-only')
% 
% figure(5);
% plot(1:48, squeeze(W5(3, 1, :))', 'b-', 1:48, squeeze(W5(3, 2, :))', 'r--', 1:48, squeeze(W4(2, 1, :))', 'k:')
% title('Forecast Error Variance Decomposition: Hours Worked')
% legend('V', 'A', 'A-only')
% 
% % output: theta22 - theta23 - theta23 + theta33 ?
% 
% %% 