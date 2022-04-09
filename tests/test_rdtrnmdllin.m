rdtrnmdllin = ss(A,B,C,D);

% Time
t_start = 0;
t_end   = 30;
t_s     = 0.1;
t_span  = (t_start:t_s:t_end)';

% Control Signals
theta = 5*(pi/180)*ones(size(t_span));

% Disturbances
dalpha_1 = zeros(size(t_span));
dalpha_2 = zeros(size(t_span));
dalpha_3 = zeros(size(t_span));
v_w2     = zeros(size(t_span));

u = [theta,dalpha_1,dalpha_2,dalpha_3,v_w2];

% Initial Conditions
v_y0   = 0;
w_zt0  = 0;
w_zs0  = 0;
Y_0    = 0;
phi_t0 = 0;
phi_s0 = 0;
x_0    = [v_y0,w_zt0,w_zs0,Y_0,phi_t0,phi_s0];

% Simulate
[~,~,x] = lsim(rdtrnmdllin,u,t_span,x_0);

v_y = x(:,1);
Y   = x(:,4);

figure;
subplot(2,1,1);
plot(t_span,v_y);
subplot(2,1,2);
plot(t_span,Y);