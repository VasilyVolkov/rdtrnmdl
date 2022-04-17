%% Prepare test scenario to compare the linear and nonlinear road train models

% Load and prepare input signals
load('rdtrnmdlcmptest.mat');

a_time = test_scenario.getElement("a").Values.Time;
a_data = test_scenario.getElement("a").Values.Data;
a      = interp1(a_time,a_data,t_span,...
    'linear','extrap');

theta_time = test_scenario.getElement("theta").Values.Time;
theta_data = test_scenario.getElement("theta").Values.Data;
theta      = interp1(theta_time,theta_data,t_span,...
    'linear','extrap');

dalpha_1_time = test_scenario.getElement("dalpha_1").Values.Time;
dalpha_1_data = test_scenario.getElement("dalpha_1").Values.Data;
dalpha_1      = interp1(dalpha_1_time,dalpha_1_data,t_span,...
    'linear','extrap');

dalpha_2_time = test_scenario.getElement("dalpha_2").Values.Time;
dalpha_2_data = test_scenario.getElement("dalpha_2").Values.Data;
dalpha_2      = interp1(dalpha_2_time,dalpha_2_data,t_span,...
    'linear','extrap');

dalpha_3_time = test_scenario.getElement("dalpha_3").Values.Time;
dalpha_3_data = test_scenario.getElement("dalpha_3").Values.Data;
dalpha_3      = interp1(dalpha_3_time,dalpha_3_data,t_span,...
    'linear','extrap');

v_w2_time = test_scenario.getElement("v_w2").Values.Time;
v_w2_data = test_scenario.getElement("v_w2").Values.Data;
v_w2      = interp1(v_w2_time,v_w2_data,t_span,...
    'linear','extrap');

u_lin  = [theta,dalpha_1,dalpha_2,dalpha_3,v_w2];
u_nlin = [a,theta,dalpha_1,dalpha_2,dalpha_3,v_w2];

% Set the initial conditions
v_x0   = 54/3.6;
v_y0   = 0;
w_zt0  = 0;
w_zs0  = 0;
X_0    = 0;
Y_0    = 0;
phi_t0 = 0;
phi_s0 = 0;

x_lin0  = [v_y0,w_zt0,w_zs0,Y_0,phi_t0,phi_s0];
x_nlin0 = [v_x0,v_y0,w_zt0,w_zs0,X_0,Y_0,phi_t0,phi_s0];