%% Road Train Parameters

% Road Train Longitudinal Speed in [m/s]
v_x = 20;

% Rolling Resistance Coefficient [null]
p.f_0 = 0.0038;

% Ambient Air Density in [kg/m3]
p.rho = 1.2;

% Air Drag Coefficient (Lateral Cross-Section) [null]
p.C_x = 0.5;
p.C_y = 0.9;

% Frontal and Lateral Cross-Sectional Areas of the Road Train in [m^2]
p.A_x  = 7.5;
p.A_y1 = 7.5;
p.A_y2 = 41.6;

% Road Train Length, Width and Height of CG in [m]
p.b = 2.045;
%l = 11.695;
%h = 1.615;

% Distances from the Truck Fifth Wheel to the Truck and Semitrailer in [m]
% CG
p.d_t = 1.815;
p.d_s = 4.34;

% Mass of the Truck and Semitrailer and Road Train in [kg]
p.m_t = 7000;
p.m_s = 15000;
p.m_rt = p.m_t + p.m_s;

% Rotating Mass factor [null]
p.sigma = 1.05;

% Moment of Inertia of the Truck and Semitrailer in [kg*m2] 
p.J_t = 30000;
p.J_s = 1/12*m_s*((2*d_s)^2 + b^2);
% J_t = 15000;
% J_s = 20000;

% Distances in [m] from the Truck and Semitrailer CGs to its axles
p.l_1 = 1.2; 
p.l_2 = 2.38;
p.l_3 = 4.34;

% Cornering Stiffness in [N/rad] Coefficients for each Road Train Axle 
p.k_1 = 50000;
p.k_2 = 150000;
p.k_3 = 3*150000;

% Gravitational Acceleration in [m/s^2]
p.g = 9.81;

% Road Train Weight, Fifth Wheel and Axle Load in [N]
P_t  = m_t*g;
P_s  = m_s*g;
P_fw = P_s*l_3/(d_s + l_3);
p.P_3  = P_s - P_fw;
p.P_2  = (P_t*l_1 + P_fw*(l_1 + d_t))/(l_1 + l_2);
p.P_1  = P_t + P_fw - P_2;