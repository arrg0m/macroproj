% flags and parameters
useunitvec = 1;
lambdaIV = 10^(-16);
diff_h = 0;
diff_r = 1;
diff_pi = 1;

% import
[tmp1, tmp2, raw] = xlsread('ishock_Data_KR.xlsx','Export','J2:O63');
R = cellfun(@(x) ~isnumeric(x) || isnan(x),raw); % Find non-numeric cells
raw(R) = {0.0}; % Replace non-numeric cells
data = cell2mat(raw);
dp = data(:,1);
da = data(:,2);
h = data(:,3);
r = data(:,4);
pi = data(:,5);
re = data(:,6); 
clearvars data raw R columnIndices;

% data sampling
year_begin = 2002;
quat_begin = 1;
begin_smpl = [2002, 2]; end_smpl = [2016, 4];

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
    
    if year == begin_smpl(1) && quat == begin_smpl(2)
        begin_itr = itr;
    end
    if year == end_smpl(1) && quat == end_smpl(2)
        end_itr = itr;
    end
end

smpl = begin_itr:end_itr;

if diff_h == 1
    h = h - [0; h(1:end-1)];
end
if diff_r == 1
    r = r - [0; r(1:end-1)];
end
if diff_pi == 1
    pi = pi - [0; pi(1:end-1)];
end

is_diff5 = [1 1 diff_h diff_r diff_pi];
is_diff4 = [1 diff_h diff_r diff_pi];

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

pi1 = [0; pi(1:end-1)];
pi2 = [0; 0; pi(1:end-2)];
pi3 = [0; 0; 0; pi(1:end-3)];
pi4 = [0; 0; 0; 0; pi(1:end-4)];
dpi = pi - pi1;
dpi1 = pi1 - pi2;
dpi2 = pi2 - pi3;
dpi3 = pi3 - pi4;

if useunitvec == 1
    unitvec = ones(len, 1);
end
%% BLOCK02 INSTRUMENTAL VARIABLE ESTIMATION
Z0 = [dp1, dp2, dp3, dp4, da1, da2, da3, da4, h1, h2, h3, h4, r1, r2, r3, r4, pi1, pi2, pi3, pi4];

% Equation 1: on dp
y = dp; y = y(smpl);
X = [dp1, dp2, dp3, dp4, dda, dda1, dda2, dda3, dh, dh1, dh2, dh3, dr, dr1, dr2, dr3, dpi, dpi1, dpi2, dpi3];
if useunitvec == 1
    X = [X, unitvec];
end
X = X(smpl,:);
Z = Z0(smpl,:);
[coef1, resid1] = myIV(y, X, Z, lambdaIV);

% Equation 2: on da
y = da; y = y(smpl);
X = [dp, dp1, dp2, dp3, dp4, da1, da2, da3, da4, dh, dh1, dh2, dh3, dr, dr1, dr2, dr3, dpi, dpi1, dpi2, dpi3]; 
if useunitvec == 1
    X = [X, unitvec];
end
X = X(smpl,:);
Z = [Z0(smpl, :), resid1];
[coef2, resid2] = myIV(y, X, Z, lambdaIV);

% Equation 3: on h
y = h; y = y(smpl);
X = [dp, dp1, dp2, dp3, da, da1, da2, da3, h1, h2, h3, r, r1, r2, r3, pi, pi1, pi2, pi3]; 
if useunitvec == 1
    X = [X, unitvec];
end
X = X(smpl,:);
Z = [Z0(smpl,:), resid1, resid2]; 
[coef3, resid3] = myIV(y, X, Z, lambdaIV);

% Equation 4: on r
y = r; y = y(smpl);
X = [dp, dp1, dp2, dp3, da, da1, da2, da3, h, h1, h2, h3, r1, r2, r3, pi, pi1, pi2, pi3];
if useunitvec == 1
    X = [X, unitvec];
end
X = X(smpl,:);
Z = [Z0(smpl,:), resid1, resid2]; 
[coef4, resid4] = myIV(y, X, Z, lambdaIV); 

% Equation 5: on pi
y = pi; y = y(smpl);
X = [dp, dp1, dp2, dp3, da, da1, da2, da3, h, h1, h2, h3, r, r1, r2, r3, pi1, pi2, pi3];
if useunitvec == 1
    X = [X, unitvec];
end
X = X(smpl,:);
Z = [Z0(smpl,:), resid1, resid2];
[coef5, resid5] = myIV(y, X, Z, lambdaIV);

% Coefficient Recovery
A5 = [1, -coef1(5), -coef1(9), -coef1(13), -coef1(17);
     -coef2(1), 1, -coef2(10), -coef2(14), -coef2(18);
     -coef3(1), -coef3(5), 1, -coef3(12), -coef3(16);
     -coef4(1), -coef4(5), -coef4(9), 1, -coef4(16);
     -coef5(1), -coef5(5), -coef5(9), -coef5(13), 1];

G51 = [coef1(1), -coef1(5) + coef1(6), -coef1(9) + coef1(10), -coef1(13) + coef1(14), -coef1(17) + coef1(18);
        coef2(2), coef2(6), -coef2(10) + coef2(11), -coef2(14) + coef2(15), -coef2(18) + coef2(19);
        -coef3(1) + coef3(2), -coef3(5) + coef3(6),  1 + coef3(9), -coef3(12) + coef3(13), -coef3(16) + coef3(17);
        -coef4(1) + coef4(2), -coef4(5) + coef4(6), -coef4(9) + coef4(10), 1 + coef4(13), -coef4(16) + coef4(17);
        -coef5(1) + coef5(2), -coef5(5) + coef5(6), -coef5(9) + coef5(10), -coef5(13) + coef5(14), 1 + coef5(17)];
    
G52 = [coef1(2), -coef1(6) + coef1(7), -coef1(10) + coef1(11), -coef1(14) + coef1(15), -coef1(18) + coef1(19);
        coef2(3), coef2(7), -coef2(11) + coef2(12), -coef2(15) + coef2(16), -coef2(19) + coef2(20);
        -coef3(2) + coef3(3), -coef3(6) + coef3(7), - coef3(9) + coef3(10), -coef3(13) + coef3(14), -coef3(17) + coef3(18);
        -coef4(2) + coef4(3), -coef4(6) + coef4(7), -coef4(10) + coef4(11), -coef4(13) + coef4(14), -coef4(17) + coef4(18);
        -coef5(2) + coef5(3), -coef5(6) + coef5(7), -coef5(10) + coef5(11), -coef5(14) + coef5(15), -coef5(17) + coef5(18)];

G53 = [coef1(3), -coef1(7) + coef1(8), -coef1(11) + coef1(12), -coef1(15) + coef1(16), -coef1(19) + coef1(20);
        coef2(4), coef2(8), -coef2(12) + coef2(13), -coef2(16) + coef2(17), -coef2(20) + coef2(21);
        -coef3(3) + coef3(4), -coef3(7) + coef3(8), -coef3(10) + coef3(11), -coef3(14) + coef3(15), -coef3(18) + coef3(19);
        -coef4(3) + coef4(4), -coef4(7) + coef4(8), -coef4(11) + coef4(12), -coef4(14) + coef4(15), -coef4(18) + coef4(19);
        -coef5(3) + coef5(4), -coef5(7) + coef5(8), -coef5(11) + coef5(12), -coef5(15) + coef5(16), -coef5(18) + coef5(19)];

G54 = [coef1(4), -coef1(8), -coef1(12), -coef1(16), -coef1(20);
        coef2(5), coef2(9), -coef2(13), -coef2(17), -coef2(21);
        -coef3(4), -coef3(8), -coef3(11), -coef3(15), -coef3(19);
        -coef4(4), -coef4(8), -coef4(12), -coef4(15), -coef4(19);
        -coef5(4), -coef5(8), -coef5(12), -coef5(16), -coef5(19)];

B5sq = cov(resid1, resid2); B5 = B5sq^0.5;
Sig5 = cov([resid1, resid2, resid3, resid4, resid5]);

cn5 = [0; 0; 0; 0; 0];
if useunitvec == 1
    cn5 = [coef1(21); coef2(22); coef3(20); coef4(20); coef5(20)];
end

%% BLOCK03 THEORETICAL IRF COMPUTATION
C5 = zeros(5, 5, 48); % coordinates: to, from, period
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
Z0 = [da1, da2, da3, da4, h1, h2, h3, h4, r1, r2, r3, r4, pi1, pi2, pi3, pi4];

% Equation 2: on da
y = da; y = y(smpl);
X = [da1, da2, da3, da4, dh, dh1, dh2, dh3, dr, dr1, dr2, dr3, dpi, dpi1, dpi2, dpi3];
if useunitvec == 1
    X = [X, unitvec];
end
X = X(smpl,:);
Z = Z0(smpl, :);
[coef2, resid2] = myIV(y, X, Z, lambdaIV);

% Equation 3: on h
y = h; y = y(smpl);
X = [da, da1, da2, da3, h1, h2, h3, r, r1, r2, r3, pi, pi1, pi2, pi3];
if useunitvec == 1
    X = [X, unitvec];
end
X = X(smpl,:);
Z = [Z0(smpl,:), resid2]; 
[coef3, resid3] = myIV(y, X, Z, lambdaIV);

% Equation 4: on r
y = r; y = y(smpl);
X = [da, da1, da2, da3, h, h1, h2, h3, r1, r2, r3, pi, pi1, pi2, pi3];
if useunitvec == 1
    X = [X, unitvec];
end
X = X(smpl,:);
Z = [Z0(smpl,:), resid2]; 
[coef4, resid4] = myIV(y, X, Z, lambdaIV); 

% Equation 5: on pi
y = pi; y = y(smpl);
X = [da, da1, da2, da3, h, h1, h2, h3, r, r1, r2, r3, pi1, pi2, pi3];
if useunitvec == 1
    X = [X, unitvec];
end
X = X(smpl,:);
Z = [Z0(smpl,:), resid2];
[coef5, resid5] = myIV(y, X, Z, lambdaIV);

A4 = [1, -coef2(5), -coef2(9), -coef2(13); 
    -coef3(1), 1, -coef3(8), -coef3(12);
    -coef4(1), -coef4(5), 1, -coef4(12);
    -coef5(1), -coef5(5), -coef5(9), 1];

G41 = [coef2(1), -coef2(5) + coef2(6), -coef2(9) + coef2(10), -coef2(13) + coef2(14); 
    -coef3(1) + coef3(2), 1 + coef3(5), -coef3(8) + coef3(9), -coef3(12) + coef3(13);
    -coef4(1) + coef4(2), -coef4(5) + coef4(6), 1 + coef4(9), -coef4(12) + coef4(13);
   -coef5(1) + coef5(2), -coef5(5) + coef5(6), -coef5(9) + coef5(10), 1 + coef5(13)];

G42 = [coef2(2), -coef2(6) + coef2(7), -coef2(10) + coef2(11), -coef2(14) + coef2(15); 
    -coef3(2) + coef3(3), -coef3(5) + coef3(6), -coef3(9) + coef3(10), -coef3(13) + coef3(14);
    -coef4(2) + coef4(3), -coef4(6) + coef4(7), -coef4(9) + coef4(10), -coef4(13) + coef4(14);
    -coef5(2) + coef5(3), -coef5(6) + coef5(7), -coef5(10) + coef5(11), -coef5(13)  + coef5(14)];

G43 = [coef2(3), -coef2(7) + coef2(8), -coef2(11) + coef2(12), -coef2(15) + coef2(16); 
    -coef3(3) + coef3(4), -coef3(6) + coef3(7), -coef3(10) + coef3(11), -coef3(14) + coef3(15);
    -coef4(3) + coef4(4), -coef4(7) + coef4(8), -coef4(10) + coef4(11), -coef4(14) + coef4(15);
    -coef5(3) + coef5(4), -coef5(7) + coef5(8), -coef5(11) + coef5(12), -coef5(14) + coef5(15)];

G44 = [coef2(4), -coef2(8), -coef2(12), -coef2(16); 
    -coef3(4), -coef3(7), -coef3(11), -coef3(15);
    -coef4(4), -coef4(8), -coef4(11), -coef4(15);
    -coef5(4), -coef5(8), -coef5(12), -coef5(15)];

B4sq = var(resid2); B4 = B4sq^0.5;
Sig4 = cov([resid2, resid3, resid4, resid5]);

cn4 = [0; 0; 0; 0];
if useunitvec == 1
    cn4 = [coef2(17); coef3(16); coef4(16); coef5(16)];
end
%% BLOCK03' THEORETICAL IRF COMPUTATION
C4 = zeros(4, 4, 48);
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

%% BLOCK04 SIMULATED IRF COMPUTATION
% http://www.dynare.org/DynareWiki/IrFs

num_sim = 100;
len_pre = 100;
len_irf = 50;
len_sim = len_pre + len_irf;
smpl_irf = (len_pre+1):len_sim;
smpl_irf_plt = 1:len_irf;

% 5-var

y0 = [0, 0, 0, 0, 0]';
y = zeros(5, len_sim, num_sim);
y1 = zeros(5, len_sim, num_sim);
y2 = zeros(5, len_sim, num_sim);

for itr = 1:num_sim
    e = Sig5^0.5 * randn(5, len_sim);
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

z0 = [0, 0, 0, 0]';
z = zeros(4, len_sim, num_sim);
z1 = zeros(4, len_sim, num_sim);

for itr = 1:num_sim
    f = Sig4^0.5 * randn(4, len_sim);
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

myirf = zeros(3, 6, len_irf);
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
myirf(1, 4, :) = y1_res(4,smpl_irf) - y_res(4,smpl_irf); % labor hour response
myirf(2, 4, :) = y2_res(4,smpl_irf) - y_res(4,smpl_irf);
myirf(3, 4, :) = z1_res(3,smpl_irf) - z_res(3,smpl_irf);
if diff_r == 1
    myirf(1, 4, :) = cumsum(myirf(1, 4, :));
    myirf(2, 4, :) = cumsum(myirf(2, 4, :));
    myirf(3, 4, :) = cumsum(myirf(3, 4, :));
end
myirf(1, 5, :) = y1_res(5,smpl_irf) - y_res(5,smpl_irf); % labor hour response
myirf(2, 5, :) = y2_res(5,smpl_irf) - y_res(5,smpl_irf);
myirf(3, 5, :) = z1_res(4,smpl_irf) - z_res(4,smpl_irf);
if diff_pi == 1
    myirf(1, 5, :) = cumsum(myirf(1, 5, :));
    myirf(2, 5, :) = cumsum(myirf(2, 5, :));
    myirf(3, 5, :) = cumsum(myirf(3, 5, :));
end
myirf(1, 6, :) = myirf(1, 2, :) + myirf(1, 3, :); % output = productivity * hours
myirf(2, 6, :) = myirf(2, 2, :) + myirf(2, 3, :);
myirf(3, 6, :) = myirf(3, 2, :) + myirf(3, 3, :);
myirf = myirf*100; % convert to percentage unit

%% BLOCK05 PLOT'em all
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
plot(smpl_irf_plt, squeeze(myirf(1,6,:)), 'b-', smpl_irf_plt, squeeze(myirf(2,6,:)), 'r--', smpl_irf_plt, squeeze(myirf(3,6,:)), 'k:');
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