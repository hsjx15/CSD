% Control System Design Project
% Authur: hshi17 10/10/18
close all; clear; clc;

deg2rad = pi/180;
rad2deg = 180/pi;

% % harsh start
% theta_init = pi;
% psi_init = pi/36;
% phi_init = pi/2;
% % observer_init = [pi/2; pi/2; pi/18; pi/2; pi; pi/2];
% % observer_init = [pi/2; 0; pi/18; 0; pi/3; 0];
% observer_init = [0; 0; 0; 0; 0; 0];

% normal start
theta_init = pi/18;
psi_init = pi/36;
phi_init = pi/18;
observer_init = [0; 0; 0; 0; 0; 0];


%% define parameters
M = 0.6;        % kg, mass of whole robot
m = 0.03;       % kg, masss of wheels
L = 0.125;      % m, distance to center of mass
R = 0.028;      % m, wheel radius
W = 0.126;      % m, body width
D = 0.07;       % m, body depth
g = 9.81;       % m/s^2, gravity acceleration
R_m = 6.83;     % Omega, resistance
J_psi = M * L^2 / 3;        % kg*m^2, moment of inertia of robot
J_phi = M * (W^2 + D^2) / 12;       % kg*m^2, moment of inertia of body yaw
J_w = m * R^2 /2;       % kg*m^2, moment of inertia of wheels
J_m = 10^(-5);      % kg*m^2, moment of inertia of DC motors
k_t = 0.305;        % Nm/A, torque constant
k_b = 0.556;        % V/(rad/s), EMF constant
f_m = 0.0022;       % friction coefficient between DC motor and body
f_w = 0;        % friction coefficient between wheel and floor
n = 1;      % gear ratio

alpha = n * k_t / R_m;
beta = n * k_t * k_b /R_m + f_m;

%% linearized state-space function coefficient matrices
m_1 = (2 * m + M) * R^2 + 2 * J_w + 2 * n^2 * J_m;
m_2 = M * L * R - 2 * n^2 * J_m;
m_3 = M * L^2 + J_psi + 2 * n^2 * J_m;
m_4 = m * W^2 /2 + J_phi + W^2 / 2 / R^2 * (J_w + n^2 * J_m);

a_22 = - (2 * beta * (m_2 + m_3) + 2 * m_3 * f_w) / (m_1 * m_3 - (m_2)^2);
a_23 = - (m_2 * M * g * L) / (m_1 * m_3 - (m_2)^2);
a_24 = 2 * beta * (m_2 + m_3) / (m_1 * m_3 - (m_2)^2);
a_42 = (2 * beta * (m_1 + m_2) + 2 * m_2 * f_w) / (m_1 * m_3 - (m_2)^2);
a_43 = m_1 * M * g * L / (m_1 * m_3 - (m_2)^2);
a_44 = - 2 * beta * (m_1 + m_2) / (m_1 * m_3 - (m_2)^2);
a_66 = - W^2 * (beta - f_w) / 2 / m_4 / R^2;

A = [0      1       0       0       0       0;
     0      a_22    a_23    a_24    0       0;
     0      0       0       1       0       0;
     0      a_42    a_43    a_44    0       0;
     0      0       0       0       0       1;
     0      0       0       0       0       a_66];
 
b_21 = alpha * (m_2 + m_3) / (m_1 * m_3 - (m_2)^2);
b_22 = alpha * (m_2 + m_3) / (m_1 * m_3 - (m_2)^2);
b_41 = - alpha * (m_1 + m_2) / (m_1 * m_3 - (m_2)^2);
b_42 = - alpha * (m_1 + m_2) / (m_1 * m_3 - (m_2)^2);
b_61 = - W * alpha / 2 / m_4 / R;
b_62 = W * alpha / 2 / m_4 / R;
 
B = [0     0;
     b_21  b_22;
     0     0;
     b_41  b_42;
     0     0;
     b_61  b_62];
  
C = [1 0 0 0 0 0;
     0 0 1 0 0 0;
     0 0 0 0 1 0];

D = zeros(3,2);

sys= ss(A, B, C, D);
%% design controller

A_con = A;
B_con = B;

if ~(length(A) - rank(ctrb(A_con,B_con)))
    disp('system is controllable');
else
    disp('system is not controllable');
end

% K = place(A, B, [-1e1, -1e1, -1, -2, - 1e3, -1e3])

% K for harsh start
% K = place(A, B, 1e2*[-3.4577 + 0.0000i, -0.0857 + 0.0000i, ...
%     -0.0757 + 0.0092i, -0.0757 - 0.0092i, ...
%     -1.8270 + 0.0000i, -0.0519 + 0.0000i])
% K = place(A, B, 1e2*[-3.5, -0.09, ...
%     -0.08 + 0.01i, -0.08 - 0.01i, ...
%     -1.83, -0.05])
% K = place(A, B, 1e2*[-3.5, -0.1, ...
%     -0.08, -0.08, ...
%     -2, -0.05])
% K = place(A, B, 1e2*[-3.5, -0.1, ...
%     -0.08, -0.08, ...
%     -2, -0.01])
% K = place(A, B, 1e2*[-1e2, -1e1, ...
%      -0.1, -8e-2, -8e-2, -5e-3])
% K = place(A, B, 1e2*[-5, -3, ...
%      -0.1, -8e-2, -8e-2, -5e-3])  % optimum psi and phi

% K for normal start
% K = place(A, B, 1e2*[-3.4577 + 0.0000i, -0.0857 + 0.0000i, ...
%     -0.0757 + 0.0092i, -0.0757 - 0.0092i, ...
%     -1.8270 + 0.0000i, -0.0519 + 0.0000i])
% K = place(A, B, 1e2*[-3.5, -0.09, ...
%     -0.08, -0.08, ...
%     -2, -5e-2])
% K = place(A, B, 1e2*[-5, -1e-2, ...
%     -1.5, -1.5, ...
%     -2, -1e-2])
% K = place(A, B, 1e2*[-5, -2e-2, ...
%     -0.8, -0.8, ...
%     -1, -2e-2])
K = place(A, B, 1e2*[-5, ...
    -1, -0.8, -0.8, ...
    -2e-2, -2e-2])

     

% alpha_con = 1/1e-2;
% [P, Lambda, K] = care(A, B*sqrt(alpha_con), C.'*C);
% K = K*sqrt(alpha_con)

% Q_con = C.'*C;
% Q_con = [
%     1e-3   0   0   0   0   0;
%     0   0   0   0   0   0;
%     0   0   1e3   0   0   0;
%     0   0   0   0   0   0;
%     0   0   0   0   1   0;
%     0   0   0   0   0   0];
% Q_con = [
%     1   0   0   0   0   0;
%     0   0   0   0   0   0;
%     0   0   1e9   0   0   0;
%     0   0   0   0   0   0;
%     0   0   0   0   1   0;
%     0   0   0   0   0   0];
% R_con = 1e-3*eye(2);
% S_con = zeros(6,2);
% E_con = eye(6);
% [P, Lambda, K] = care(A, B, Q_con, R_con, S_con, E_con)

% Q_lqr = [
%     1e-6   0   0   0   0   0;
%     0   0   0   0   0   0;
%     0   0   1e9   0   0   0;
%     0   0   0   0   0   0;
%     0   0   0   0   1   0;
%     0   0   0   0   0   0];
% R_lqr = [
%     1e-3   0;
%     0   1e-3];
% N_lqr = zeros(6, 2)
% [K, S_con, eval_con] = lqr(A, B, Q_lqr, R_lqr, N_lqr)

%% design observer

C_obsv = eye(6);
D_obsv = zeros(6, size(B,2)+size(C_obsv,1));

if ~(length(A) - rank(ctrb(A.', C_obsv.')))
    disp('system is observable');
else
    disp('system is not observable');
end

% L_obsv = place(A', C_obsv', [-5, -10, -1e2, -1e2, -1e4, -1e4])

alpha_obsv = 1/1e-2;
[Q, Lambda_obsv, L_obsv] = care(A.', C_obsv.'*sqrt(alpha_obsv), B*B.');
L_obsv = L_obsv * alpha_obsv

%% discrete system
sys = ss(A, B, C, D);
h_sample = 4e-3;
sysd = c2d(sys, h_sample);

Q_lqrd = eye(6);
Q_lqrd(2,2) = 0.5;
Q_lqrd(3,3) = 100;
R_lqrd = eye(2);
R_lqrd = R_lqrd/15;
[K_con_dis,~,~] = lqrd(A, B, Q_lqrd, R_lqrd, 4e-3)

% alpha_con_dis = 1e-1;
% [P_con_dis, Lambda_con_dis, K_con_dis] = dare(sysd.A, sysd.B*sqrt(alpha_con_dis), sysd.C.'*sysd.C);
% K_con_dis = K_con_dis*sqrt(alpha_con_dis)

%% disgonosis
% when matrix M for (theta, psi) is not invertible
% cos_phi = (sqrt(m_1 * m_3) + 2 * n^2 * J_m) / (M * L * R)

% Lin_Segway = linmod('CSD_Segway_Nonlinear');
% Lin_Segway.a
% A