% Time
t_start = 0;
t_end   = 30;
t_s     = 0.1;
t_span  = (t_start:t_s:t_end)';

% Control Signals
a     = zeros(size(t_span));
theta = 5*(pi/180)*ones(size(t_span));

% Disturbances
dalpha_1 = zeros(size(t_span));
dalpha_2 = zeros(size(t_span));
dalpha_3 = zeros(size(t_span));
v_w2     = zeros(size(t_span));

u        = [a,theta,dalpha_1,dalpha_2,dalpha_3,v_w2];

% Initial Conditions
v_x0   = 54/3.6;
v_y0   = 0;
w_zt0  = 0;
w_zs0  = 0;
X_0    = 0;
Y_0    = 0;
phi_t0 = 0;
phi_s0 = 0;
x_0    = [v_x0,v_y0,w_zt0,w_zs0,X_0,Y_0,phi_t0,phi_s0];

opts  = odeset(...
    'RelTol',1e-2,...
    'AbsTol',1e-4);
[t,x] = ode45(@(t,x) rdtrnmdlnlin(t,x,t_span,u,p),t_span,x_0,opts);

v_x   = x(:,1);
v_y   = x(:,2);
w_zt  = x(:,3);
w_zs  = x(:,4);
X     = x(:,5);
Y     = x(:,6);
phi_t = x(:,7);
phi_s = x(:,8);

figure
subplot(2,1,1)
plot(t,v_x);
subplot(2,1,2)
plot(t,v_y);

figure
plot(X,Y)