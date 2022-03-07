%% Road Train Parameters

% Road Train Longitudinal Speed in [m/s]
v_x = 20;

% Rolling Resistance Coefficient [null]
f_0 = 0.0038;

% Ambient Air Density in [kg/m3]
rho = 1.2;

% Air Drag Coefficient (Lateral Cross-Section) [null]
C_x = 0.5;
C_y = 0.9;

% Frontal and Lateral Cross-Sectional Areas of the Road Train in [m^2]
A_x   = 7.5;
A_y_1 = 7.5;
A_y_2 = 41.6;

% Road Train Length, Width and Height of CG in [m]
b = 2.045;
%l = 11.695;
%h = 1.615;

% Distances from the Truck Fifth Wheel to the Truck and Semitrailer in [m]
% CG
d_1 = 1.815;
d_2 = 4.34;

% Mass of the Truck and Semitrailer in [kg]
m_t = 7000;
m_s = 15000;

% Rotating Mass factor [null]
sigma = 1.05;

% Moment of Inertia of the Truck and Semitrailer in [kg*m2] 
J_t = m_t*d_1^2 + 30000;
J_s = m_s*d_2^2 + 1/12*m_s*((2*d_2)^2 + b^2);
% J_t = 15000;
% J_s = 20000;

% Distances in [m] from the Truck and Semitrailer CGs to its axles
l_1 = 1.2; 
l_2 = 2.38;
l_3 = 4.34;


% Cornering Stiffness in [N/rad] Coefficients for each Road Train Axle 
k_1 = 50000;
k_2 = 150000;
k_3 = 3*150000;

% Gravitational Acceleration in [m/s^2]
g = 9.81;

% Road Train Weight in [N]
P_t  = m_t*g;
P_s  = m_s*g;

% Road Train Fifth Wheel Load in [N]
P_fw = P_s*l_3/(d_2 + l_3);

% Road Train Axle Load in [N]
P_3 = P_s - P_fw;
P_2 = (P_t*l_1 + P_fw*(l_1 + d_1))/(l_1 + l_2);
P_1 = P_t + P_fw - P_2;

