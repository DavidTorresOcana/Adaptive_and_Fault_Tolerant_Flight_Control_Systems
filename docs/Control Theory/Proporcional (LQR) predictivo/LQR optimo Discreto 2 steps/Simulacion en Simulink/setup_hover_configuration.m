% SETUP_HOVER_CONFIGURATION
%
% SETUP_HOVER_CONFIGURATION sets and returns the model model parameters 
% of the Quanser 3 DOF Hover plant.
%
%
% Copyright (C) 2010 Quanser Consulting Inc.
% Quanser Consulting Inc.
%
%
function [ Ktn, Ktc, Kf, l, Jy, Jp, Jr, g ] = setup_hover_configuration( )
%
% Gravitational Constant (m/s^2)
g = 9.81;
% Motor Armature Resistance (Ohm)
Rm = 0.83;
% Motor Current-Torque Constant (N.m/A)
Kt_m = 0.0182;
% Motor Rotor Moment of Inertia (kg.m^2)
Jm = 1.91e-6;
% Moving Mass of the Hover system (kg)
m_hover = 2.85;
% Mass of each Propeller Section = motor + shield + propeller + body (kg)
m_prop = m_hover / 4;
% Distance between Pivot to each Motor (m)
l = 7.75*0.0254;
% Propeller Force-Thrust Constant found Experimentally (N/V)
Kf = 0.1188;
% Propeller Torque-Thrust Constant found Experimentally (N.m/V)
Kt_prop = 0.0036;
% Normal Rotation Propeller Torque-Thrust Constant  (N.m/V)
Ktn = Kt_prop;
% Counter Rotation Propeller Torque-Thrust Constant  (N.m/V)
Ktc = -Kt_prop;
% Equivalent Moment of Inertia of each Propeller Section (kg.m^2)
Jeq_prop = Jm + m_prop*l^2;
% Equivalent Moment of Inertia about each Axis (kg.m^2)
Jp = 2*Jeq_prop;
Jy = 4*Jeq_prop;
Jr = 2*Jeq_prop;
%
% end of setup_hover_configuration()