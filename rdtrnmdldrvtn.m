%% Info
 
% This script derives the equations of motion for the road train planar
% motion with semitrailer. The model has 8 states and 2 inputs. The
% equations are obtained with the Symbolic Math Toolbox.
%
% Copyright 2022 Vasily Volkov, Gulshat Galimova
 
%% Cleanup

clear;
clc;

%% Declare symbolic variables
 
% Model parameters
syms d_t d_s   l_1  l_2   l_3 k_1 k_2 k_3;
syms m_t m_s   J_t  J_s sigma   g f_0;
syms C_x C_y A_xrt A_yt  A_ys rho;
p = [d_t d_s   l_1  l_2   l_3 k_1 k_2 k_3 ...
     m_t m_s   J_t  J_s sigma   g f_0 ...
     C_x C_y A_xrt A_yt  A_ys rho].';
assume(in(p,'real') & (p > 0));
 
% Model states
syms v_x v_y w_zt w_zs X Y phi_t phi_s;
x = [v_x,v_y,w_zt,w_zs,X,Y,phi_t,phi_s].';
assume(in(x,'real'));
 
% Model control inputs
syms a theta;
u = [a,theta].';
assume(in(u,'real'));

% Wind Speed
syms v_w2;
assume(in(v_w2,'real'));

% Angles
syms dphi alpha_1 alpha_2 alpha_3 dalpha_1 dalpha_2 dalpha_3;
angles = [dphi,alpha_1,alpha_2,alpha_3,dalpha_1,dalpha_2,dalpha_3].';
assume(in(angles,'real') & (angles >= -pi) & (angles < pi));

% Total Mass
syms m_rt;
assume(in(m_rt,'real') & (m_rt > 0));

% Road Train Weight, Fifth Wheel and Axle Load
syms P_t P_s P_fw P_1 P_2 P_3;
loads = [P_t,P_s,P_fw,P_1,P_2,P_3].';
assume(in(loads,'real') & (loads > 0));

% Road Train Axle Reactive and Drag Forces
syms R_x1 R_x2 R_x3 R_y1 R_y2 R_y3 F_wxrt F_wyt F_wys;
forces = [R_x1,R_x2,R_x3,R_y1,R_y2,R_y3,F_wxrt,F_wyt,F_wys].';
assume(in(forces,'real'));

%% Derive the equations of motion
 
% Mass Matrix
M_1 = [
       sigma*m_rt        0                 0 -m_s*d_s*sin(dphi);
                0     m_rt           m_t*d_t -m_s*d_s*cos(dphi);
                0  m_t*d_t J_t + m_t*d_t*d_t                  0;
                0 -m_s*d_s                 0  J_s + m_s*d_s*d_s;
      ];
M_2 = eye(4);
M   = blkdiag(M_1,M_2);
 
% Forces
Q_1 = m_rt*v_y*w_zt + m_t*d_t*w_zt*w_zt - m_s*d_s*w_zs*w_zs*cos(dphi) - ...
    F_wxrt + ...
    R_x1*cos(theta) + R_x2 - ...
    R_y1*sin(theta) + R_x3*cos(dphi) + R_y3*sin(dphi);
Q_2 = -sigma*m_rt*v_x*w_zt + m_s*d_s*w_zs*w_zs*sin(dphi) + ...
    F_wyt + F_wys + ...    
    R_x1*sin(theta) + R_y1*cos(theta) + ...
    R_y2 - R_x3*sin(dphi) + R_y3*cos(dphi);
Q_3 = (sigma - 1)*m_rt*v_x*v_y - m_t*d_t*v_x*w_zt + ...
    F_wyt*d_t + ...
    R_x1*(d_t + l_1)*sin(theta) + R_y1*(d_t + l_1)*cos(theta) + ...
    R_y2*(d_t - l_2);
Q_4 = m_s*d_s*v_x*w_zt - ...
    F_wys*d_s + ...
    R_y3*(d_s + l_3);
Q_5 = v_x*cos(phi_t) - v_y*sin(phi_t);
Q_6 = v_x*sin(phi_t) + v_y*cos(phi_t);
Q_7 = w_zt;
Q_8 = w_zs;
Q   = simplify([Q_1;Q_2;Q_3;Q_4;Q_5;Q_6;Q_7;Q_8]);

%% Simplify formulae of nonlinear model

M_inv  = inv(M);
M_invQ = [M_inv,Q];

[M_invQ,p_0] = subexpr(collect(M_invQ),'p_0');
[M_invQ,p_1] = subexpr(collect(M_invQ),'p_1');
[M_invQ,p_2] = subexpr(collect(M_invQ),'p_2');
[M_invQ,p_3] = subexpr(collect(M_invQ),'p_3');
[M_invQ,p_4] = subexpr(collect(M_invQ),'p_4');

% Get simplified matrices
M_inv = M_invQ(1:size(M,1),1:size(M,2));
Q     = M_invQ(1:(size(Q,1)),(size(M,2) + 1):(size(M,2)+size(Q,2)));
