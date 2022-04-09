%% Road Train Linear Model

% Road Train Longitudinal Speed in [m/s]
v_x = 54/3.6;

% Inertia Matrix
M = [
      p.m_t+p.m_s	      p.m_t*p.d_t        -p.m_s*p.d_s;
      p.m_t*p.d_t p.J_t+p.m_t*p.d_t^2                   0;
     -p.m_s*p.d_s                   0 p.J_s+p.m_s*p.d_s^2;    
    ];
M = blkdiag(M,eye(3));

% State-to-State Matrix
A = zeros(6);
A(:,1) = [
                                                               -(p.k_1+p.k_2+p.k_3)/v_x;
          (p.sigma-1)*(p.m_t+p.m_s)*v_x - (p.k_1*(p.d_t+p.l_1)+p.k_2*(p.d_t-p.l_2))/v_x;
                                                              (p.k_3*(p.d_s+p.l_3))/v_x;
                                                                                      1;
                                                                                      0;
                                                                                      0;
         ];
 
A(:,2) = [
                                 -p.sigma*(p.m_rt)*v_x - (p.k_1*p.l_1-p.k_2*p.l_2)/v_x;
          -(p.k_1*p.l_1*(p.d_t+p.l_1)-p.k_2*p.l_2*(p.d_t-p.l_2))/v_x - p.m_t*p.d_t*v_x;
                                                                       p.m_s*p.d_s*v_x;
                                                                                     0;
                                                                                     1;
                                                                                     0;
         ];
     
A(:,3) = [
                         (p.k_3*p.l_3)/v_x;
                                         0;
          -(p.k_3*p.l_3*(p.d_s+p.l_3))/v_x;
                                         0;
                                         0;
                                         1;
         ];

A(:,4) = [
          0;
          0;
          0;
          0;
          0;
          0;
         ];

A(:,5) = [
          p.f_0*p.P_3 - p.k_3;
                            0;
          p.k_3*(p.d_s+p.l_3);
                          v_x;
                            0;
                            0;
         ];

A(:,6) = [
             p.k_3 - p.f_0*p.P_3;
                               0;
          -(p.k_3*(p.d_s+p.l_3));
                               0; 
                               0;
                               0;
         ];

% Control-to-State Matrix
% u = B_1*theta,
% where theta - wheel angle in [rad]
B_1 = [
                     p.k_1 - p.f_0*p.P_1; 
       (p.k_1-p.f_0*p.P_1)*(p.d_t+p.l_1);
                                       0;
                                       0;
                                       0;
                                       0;
      ];

% Disturbance-to-State Matrix
% w = B_2*[dalpha_1; dalpha_2; dalpha_3; v_w^2],
% where dalpha_i - Road Train i-th Axle Slip Angle Disturbance in [rad], 
% v_w - Side Wind Speed in [m/s]
B_2 = [
                     p.k_1,               p.k_2,                p.k_3, 0.5*p.rho*p.C_y*(p.A_yt+p.A_ys);
       p.k_1*(p.d_t+p.l_1), p.k_2*(p.d_t-p.l_2),                    0,    0.5*p.rho*p.C_y*p.A_yt*p.d_t;
                         0,                   0, -p.k_3*(p.d_s+p.l_3),   -0.5*p.rho*p.C_y*p.A_ys*p.d_s;
                         0,                   0,                    0,                               0;
                         0,                   0,                    0,                               0;
                         0,                   0,                    0,                               0;
      ];

% Input-to-State Matrix
B = [B_1, B_2];

% State-to-Output Matrix
C = [
     1 0 0 0 0 0;
     0 0 0 1 0 0;
    ];

% Input-to-Output Matrix
D = zeros(size(C,1),size(B,2));

% Linear Functional Matrix
% g = K*x = [alpha_2; alpha_3; dphi].
K = [
     -1/v_x p.l_2/v_x         0  0  0  0;
     -1/v_x         0 p.l_3/v_x  0  0  0;
          0         0         0  0  1 -1;
    ];

% Solve the state-space equations with respect to the derivatives
A = M\A;
B = M\B;
