clc;
clear;
disp('Take-home exam problem 1');

disp('Hit any key when ready...');
pause;

% Setting parameters:
chi=1;
theta=2.9869;
alfa=.661;
delta=.021;
betta=.9926;
gama=.004;
sigma_eps = .0115;
x_bar=exp(gama);
% Calculating the steady state:

R_bar=1/betta;
k_bar=8.052;
y_bar=.768;
c_bar=.567;
i_bar=.201;
N_bar=.2305;
L_bar=1-N_bar;

% Declaring the matrices. 


VARNAMES = ['capital    ',
            'consumption',
            'output     ',
            'labor      ',
            'interest   ',
            'investment ',
            'technology '];


% for k(t):
AA = [ 0
       - k_bar
       0
       0
       -(1-alfa)*(1-alfa)*y_bar/k_bar*exp((alfa-1)*gama) ];

% for k(t-1):
BB = [ 1-alfa
       (1-delta)*k_bar*(1-gama)
       0
       0
       -alfa*(1-alfa)*y_bar/k_bar*exp((alfa-1)*gama) ];

%Order:   consumption  output      labor     interest  investment
CC = [   0,            -1,         alfa,        0,      0 % Equ. 1)
          0,           0,          0,           0,        i_bar  % Equ. 2)
          -c_bar,      y_bar,      0,           0,        -i_bar      % Equ. 3)      
          chi,        -1,          1-N_bar/L_bar, 0,        0      % Equ. 4)
          0,           (1-alfa)*y_bar/k_bar*exp((alfa-1)*gama), 0,     - R_bar,  0 ];   % Equ. 5)

DD = [ alfa-1
       (1-delta)*k_bar*(1-gama)
       0
       0
       (1-alfa)*y_bar/k_bar*exp((alfa-1)*gama)*(alfa-1)-(1-delta)*exp(-gama) ];

FF = [ 0 ];

GG = [ 0 ];

HH = [ 0 ];

JJ = [ -1,  0,  0,  1,  0];

KK = [ 1,   0,  0,  0,  0];

LL = [ 0 ];

MM = [ 0 ];

NN = [0];

Sigma = [ sigma_eps^2  ];

% Setting the options:

[l_equ,m_states] = size(AA);
[l_equ,n_endog ] = size(CC);
[l_equ,k_exog  ] = size(DD);

  
PERIOD     = 4;  % number of periods per year, i.e. 12 for monthly, 4 for quarterly
GNP_INDEX  = 3; % Index of output among the variables selected for HP filter
IMP_SELECT = [1:7];
   %  a vector containing the indices of the variables to be plotted
DO_SIMUL   = 1; % Calculates simulations
SIM_LENGTH = 150;
DO_MOMENTS = 1; % Calculates moments based on frequency-domain methods
HP_SELECT  = 1:(m_states+n_endog+k_exog); % Selecting the variables for the HP Filter calcs.
% DO_COLOR_PRINT = 1;


% Starting the calculations:

do_it;


 
       

