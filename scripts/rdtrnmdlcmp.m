%% Road Train Linear and Nonlinear Model Comparison

% Set the time span of the simulation
t_s     = 0.01;
t_start = 0;
t_end   = 12;
t_span  = (t_start:t_s:t_end-t_s)';

% Prepare the comparison scenario
rdtrnmdlcmptestprep;

% Simulate the linear and nonlinear road train model responses
rdtrnmdllin     = ss(A,[B_1,B_2],C,[D_1,D_2]);
[~,t_lin,x_lin] = lsim(rdtrnmdllin,u_lin,t_span,x_lin0);
[t_nlin,x_nlin] = ...
    ode45(@(t,x) rdtrnmdlnlin(t,x,t_span,u_nlin,p),t_span,x_nlin0);

% Set the road train states under consideration
num_st    = 5;
x_lin_st  = [1,2,3,5,6];
x_nlin_st = [2,3,4,7,8];

% Select the results to plot
rdtrnmdlcmprslts;

% Plot the Results
rdtrnmdlcmprsltsplot;