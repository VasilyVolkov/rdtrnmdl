function dxdt = rdtrnmdlnlin(t,x,t_span,u,p)
%RDTRNMDLNLIN Summary of this function goes here
%   Detailed explanation goes here

%% Initialization

% Interpolate External Time-Dependent Model Inputs and Disturbances
u = interp1(t_span,u,t);

% Initialize Model States
v_x   = x(1);
v_y   = x(2);
w_zt  = x(3);
w_zs  = x(4);
X     = x(5);
Y     = x(6);
phi_t = x(7);
phi_s = x(8);

% Initialize Model Inputs
a     = u(1);
theta = u(2);

% Initialize Model Disturbances
dalpha_1 = u(3);
dalpha_2 = u(4);
dalpha_3 = u(5);
v_w2     = u(6);

% Initialize Model Parameters
f_0   = p.f_0;
rho   = p.rho;
C_x   = p.C_x;
C_y   = p.C_y;
A_xrt = p.A_xrt;
A_yt  = p.A_yt;
A_ys  = p.A_ys;
d_t   = p.d_t;
d_s   = p.d_s;
m_t   = p.m_t;
m_s   = p.m_s;
m_rt  = p.m_rt;
sigma = p.sigma;
J_t   = p.J_t;
J_s   = p.J_s;
l_1   = p.l_1; 
l_2   = p.l_2;
l_3   = p.l_3;
k_1   = p.k_1;
k_2   = p.k_2;
k_3   = p.k_3;
g     = p.g;
P_1   = p.P_1;
P_2   = p.P_2;
P_3   = p.P_3;

%% Preliminary Calculations

% Jackknifing Angle
dphi = phi_t - phi_s;

% Slip angles
alpha_1 = theta - atan2(v_y + w_zt*l_1,v_x);
alpha_2 = atan2(w_zt*l_2 - v_y,v_x);
alpha_3 = atan2(w_zs*l_3 - v_x*sin(dphi) - v_y*cos(dphi),...
               v_x*cos(dphi) - v_y*sin(dphi));
 
% Longitudinal forces
R_x1 = -f_0*P_1;
R_x2 = -f_0*P_2 + m_rt*sigma*a;
R_x3 = -f_0*P_3;
 
% Lateral forces
R_y1 = k_1*(alpha_1 + dalpha_1);
R_y2 = k_2*(alpha_2 + dalpha_2);
R_y3 = k_3*(alpha_3 + dalpha_3);
 
% Drag Force
F_wxrt = 0.5*rho*C_x*A_xrt*v_x*v_x;
F_wyt  = 0.5*rho*C_y*A_yt*v_w2;
F_wys  = 0.5*rho*C_y*A_ys*v_w2;

% Intermediate Parameters
p_0 = 1/(cos(dphi)*d_s^2*d_t^2*m_s^2*m_t + d_s^2*d_t^2*m_s*m_t^2 - ...
    m_rt*d_s^2*d_t^2*m_s*m_t + ...
    J_t*cos(dphi)*d_s^2*m_s^2 - J_t*m_rt*d_s^2*m_s + ...
    J_s*d_t^2*m_t^2 - J_s*m_rt*d_t^2*m_t - J_s*J_t*m_rt);
p_1 = 1/(sigma*d_s^2*d_t^2*m_rt^2*m_s*m_t - ...
    sigma*cos(dphi)*d_s^2*d_t^2*m_rt*m_s^2*m_t - ...
    sigma*d_s^2*d_t^2*m_rt*m_s*m_t^2 + ...
    J_t*sigma*d_s^2*m_rt^2*m_s - J_t*sigma*cos(dphi)*d_s^2*m_rt*m_s^2 + ...
    J_s*sigma*d_t^2*m_rt^2*m_t - J_s*sigma*d_t^2*m_rt*m_t^2 + ...
    J_s*J_t*sigma*m_rt^2);
p_2 = -d_t^2*m_t^2 + m_rt*d_t^2*m_t + J_t*m_rt;
p_3 = p_0*(d_t*m_s*m_t*d_s^2 + J_s*d_t*m_t);
p_4 = d_s*m_s*m_t*d_t^2 + J_t*d_s*m_s;

%% Calculation of State Derivatives

% Calculate the Inverse of the Inertia Matrix M
M_inv(1,1) = 1/(m_rt*sigma);
M_inv(1,2) = d_s*m_s*p_1*p_4*sin(dphi);
M_inv(1,3) = -d_s^2*d_t*m_s^2*m_t*p_1*sin(dphi);
M_inv(1,4) = d_s*m_s*p_1*p_2*sin(dphi);
M_inv(2,1) = 0;
M_inv(2,2) = -p_0*(m_s*m_t*d_s^2*d_t^2 + ...
    J_t*m_s*d_s^2 + J_s*m_t*d_t^2 + J_s*J_t);
M_inv(2,3) = p_3;
M_inv(2,4) = -p_0*(d_s*m_s*m_t*cos(dphi)*d_t^2 + J_t*d_s*m_s*cos(dphi));
M_inv(3,1) = 0;
M_inv(3,2) = p_3;
M_inv(3,3) = -p_0*(- cos(dphi)*d_s^2*m_s^2 + m_rt*d_s^2*m_s + J_s*m_rt);
M_inv(3,4) = d_s*d_t*m_s*m_t*p_0*cos(dphi);
M_inv(4,1) = 0;
M_inv(4,2) = -p_0*p_4;
M_inv(4,3) = d_s*d_t*m_s*m_t*p_0;
M_inv(4,4) = -p_0*p_2;
M_inv      = blkdiag(M_inv,eye(4));

% Calculate the Forces Matrix Q
Q(1,1) = d_t*m_t*w_zt^2 + m_rt*v_y*w_zt - d_s*m_s*cos(dphi)*w_zs^2 - ...
    F_wxrt + R_x1*cos(theta) + R_x2 + R_x3*cos(dphi) - ...
    R_y1*sin(theta) + R_y3*sin(dphi);
Q(2,1) = d_s*m_s*sin(dphi)*w_zs^2 + F_wyt + F_wys + ...
    R_x1*sin(theta) - R_x3*sin(dphi) + ...   
    R_y1*cos(theta) + R_y2 + R_y3*cos(dphi) - ...
    m_rt*sigma*v_x*w_zt;
Q(3,1) = F_wyt*d_t + (-d_t*m_t*v_x)*w_zt + ...
    R_x1*sin(theta)*(d_t + l_1) + R_y1*cos(theta)*(d_t + l_1) + ...
    m_rt*v_x*v_y*(sigma - 1);
Q(4,1) = d_s*m_s*v_x*w_zt - R_y3*(d_s + l_3) - F_wys*d_s;
Q(5,1) = v_x*cos(phi_t) - v_y*sin(phi_t);
Q(6,1) = v_y*cos(phi_t) + v_x*sin(phi_t);
Q(7,1) = w_zt;
Q(8,1) = w_zs;

% Calculate the Derivatives of the Model States
dxdt = M_inv*Q;

end

