% VERSION 2.0, MARCH 1997, COPYRIGHT H. UHLIG.
% EXAMPL1.M calculates through Hansens benchmark real business
% cycle model in H. Uhlig, "A toolkit solving nonlinear dynamic stochastic models easily".
% First, parameters are set and the steady state is calculated. Next, the matrices are
% declared.  In the last line, the model is solved and analyzed by calling DO_IT.M

% Copyright: H. Uhlig.  Feel free to copy, modify and use at your own risk.
% However, you are not allowed to sell this software or otherwise impinge
% on its free distribution.

% Modified by Carlos Lezama and Santiago Payró from Felipe Meza's example
% Dynamic Macroeconomics II
% ITAM - Fall 2022

% -------------------------------------------------------------------------
% Reset

clc; clear all; close all;

% -------------------------------------------------------------------------
% Parameters

N_bar = 1.0/3; % steady state employment is a third of total time endowment
Z_bar = 1; % normalization
alpha = .36; % capital share
delta = .025; % depreciation rate capital
R_bar = 1.01; % one percent real interest per quarter
psi = .95; % autocorrelation of technology shock
sigma_eps = .712; % standard deviation of technology shock (percent)
nu = 1.5; % labor supply elasticity of substitution

% -------------------------------------------------------------------------
% Steady state and other necessary variables

beta = 1 / R_bar;
YK_bar = (R_bar - 1 + delta) / alpha;
K_bar = (Z_bar / YK_bar) ^ (1 - alpha) * N_bar;
I_bar = delta * K_bar;
Y_bar = YK_bar * K_bar;
C_bar = Y_bar - I_bar;
tau = (1 - alpha) * Y_bar / (N_bar ^ nu);
L_bar = 1 / (C_bar - ((tau / nu) * (N_bar ^ nu)));
I_bar = Y_bar - C_bar;
rho = (alpha / R_bar) * (Y_bar / K_bar);
eta = 1 / (C_bar - ((tau / nu) * (N_bar ^ nu)));
theta = eta * tau * (N_bar ^ nu);
phi = (((nu * C_bar) - (tau * N_bar * nu)) ^ 2 * (1 - alpha)) / (nu ^ 3 * tau * C_bar * (N_bar ^ nu));

% -------------------------------------------------------------------------
% Matrices

VARNAMES = ['Capital    '
            'Consumption'
            'Output     '
            'Labor      '
            'Investment '
            'Technology '];

% k(t):
AA = [0
      -K_bar
      0
      0];

% k(t-1):
BB = [0
      (1 - delta) * K_bar
      alpha
      0];

% c(t) y(t) n(t) i(t)
CC = [-C_bar, Y_bar, 0, -I_bar % eq. 1
      0, 0, 0, I_bar % eq. 2
      0, -1, 1 - alpha, 0 % eq. 3
      (1 / nu) - (phi * eta * L_bar * C_bar), phi * L_bar, phi * theta * L_bar, 0]; % eq. 4

% z(t)
DD = [0
      0
      1
      0];

% k(t+1)
FF = [0];

% k(t)
GG = [-rho];

% k(t-1)
HH = [0];

% c(t+1) y(t+1) n(t+1) i(t+1)
JJ = [-eta * C_bar, rho, eta * tau * (N_bar ^ nu), 0];

% c(t) y(t) n(t) i(t)
KK = [eta * C_bar, 0, -eta * tau * (N_bar ^ nu), 0];

% z(t+1)
LL = [0];

% z(t)
MM = [0];

% AR z(t)
NN = [psi];

% AR variance
Sigma = [sigma_eps ^ 2];

% -------------------------------------------------------------------------
% Options

l_equ = size(AA, 1);
m_states = size(AA, 2);
n_endog = size(CC, 2);
k_exog = size(DD, 2);

% -------------------------------------------------------------------------
% Calculations

warnings = [];
options;
solve;
sol_out;
%do_it;

% -------------------------------------------------------------------------
% Impulse responses

T = 150 + 1;
Tend = 33;

time_count = 0;

for shock_counter = SELECT_SHOCKS % 1 : k_exog,
    response = zeros(m_states + n_endog + k_exog, HORIZON);
    response(m_states + n_endog + shock_counter, 1) = 1;
    left = [PP, zeros(m_states, n_endog), zeros(m_states, k_exog)
            RR, zeros(n_endog, n_endog), zeros(n_endog, k_exog)
            zeros(k_exog, (m_states + n_endog)), NN];
    right = eye(m_states + n_endog + k_exog) + ...
        [zeros(m_states, (m_states + n_endog)), QQ
     zeros(n_endog, (m_states + n_endog)), SS
     zeros(k_exog, (m_states + n_endog)), zeros(k_exog, k_exog)];
    response(:, 1) = right * response(:, 1);

    for time_counter = 2:T
        response(:, time_counter) = right * left * response(:, time_counter - 1);
    end

end

response = [zeros(6, 1) response];

figure(1)

for i = 1:6
    subplot(3, 2, i)
    plot(-1:T - 1, response(i, :))
    title(VARNAMES(i, :))
    xlim([-1 Tend])
end

% y(t) first 5 responses ahead
response(3, 1:6)
