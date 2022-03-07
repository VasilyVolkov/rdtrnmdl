%% Road Train Linear Model

% Inertia Matrix
M = [
      m_1+m_2	    m_1*d_1      -m_2*d_2;
      m_1*d_1 J_1+m_1*d_1^2             0;
     -m_2*d_2             0 J_2+m_2*d_2^2;    
    ];
M = blkdiag(M,eye(3));

% State-to-State Matrix
A = zeros(6);
A(:,1) = [
                                               -(k_1+k_2+k_3)/v_x;
          (s-1)*(m_1+m_2)*v_x - (k_1*(d_1+l_1)+k_2*(d_1-l_2))/v_x;
                                              (k_3*(d_2+l_3))/v_x;
                                                                1;
                                                                0;
                                                                0;
         ];
 
A(:,2) = [
                          -s*(m_1+m_2)*v_x - (k_1*l_1-k_2*l_2)/v_x;
          -(k_1*l_1*(d_1+l_1)-k_2*l_2*(d_1-l_2))/v_x - m_1*d_1*v_x;
                                                       m_2*d_2*v_x;
                                                                 0;
                                                                 1;
                                                                 0;
         ];
     
A(:,3) = [
                     (k_3*l_3)/v_x;
                                 0;
          -(k_3*l_3*(d_2+l_3))/v_x;
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
          f_0*N_3 - (k_3);
                        0;
            k_3*(d_2+l_3);
                      v_x;
                        0;
                        0;
         ];

A(:,6) = [
             k_3 - f_0*N_3;
                         0;
          -(k_3*(d_2+l_3));
                         0; 
                         0;
                         0;
         ];

% Input-to-State Matrix
B = [
                 k_1-f_0*N_1;
     (k_1-f_0*N_1)*(d_1+l_1);
                           0;
                           0;
                           0;
                           0;
    ];

% State-to-Output Matrix
C = [
     1 0 0 0 0 0;
     0 0 0 1 0 0;
    ];

% Input-to-Output Matrix
D = zeros(size(C,1),size(B,2));

% Disturbance-to-State Matrix
% w = F*[dalpha_1; dalpha_2; dalpha_3; v_w^2],
% where dalpha_i - Road Train i-th Axle Slip Angle Disturbance in [rad], 
% v_w - Side Wind Speed in [m/s]
F = [
               k_1,           k_2,            k_3, 0.5*rho*C_y*(A_y1+A_y2);
     k_1*(d_1+l_1), k_2*(d_1-l_2),              0,        0.5*rho*C_y*A_y1;
                 0,             0, -k_3*(d_2+l_3),       -0.5*rho*C_y*A_y2;
                 0,             0,              0,                       0;
                 0,             0,              0,                       0;
                 0,             0,              0,                       0;
    ];

% Linear Functional Matrix
% g = K*x = [alpha_2; alpha_3; dphi].
K = [
     -1/v_x  l_2/v_x         0  0  0  0;
     -1/v_x        0   l_3/v_x  0  0  0;
          0        0         0  0  1 -1;
    ];

% Solve the state-space equations with respect to the derivatives
A = M\A;
B = M\B;
F = M\F;
