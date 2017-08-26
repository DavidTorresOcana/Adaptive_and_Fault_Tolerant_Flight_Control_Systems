%SetupCompleto Es la funcion que define TODOS los parametros del sistema

global  v_i0 options m M Ix Iy Iz l h Kp Kb Li Ri R b c sigma a theta_0 theta_1 C_d0 A_x A_y A_z C_x C_y C_z g rho  J_T r_cg C_ax C_ay C_az
options =optimset('Display','off'); % Settings del fsolve de la velocidad inducida
%% MASAS, momentos de inercia y longitudes
% m=2.85; % Masa del cuadricoptero 2.85Kg
% m=1.67; % Equilibrado
m=1.4;
M=0.005; % Masa del disco equivalente a la helice girando. Hay que sintonizar esta masa para que el modelo Semicompleto tenga el mismo transitorio que los motores reales
Ix=0.0552;
Iy=0.0552;
Iz=0.11;
l=0.197;
h= 0.03;% Distancia de las helices al CG ( O al punto de apoyo en su caso)
h_cg=0.10; %Distancia en ejes cuerpo del punto de giro al CG
r_cg=[0;0;h_cg];
%% Parametros de los motores/helice
Kp=1.82E-2;
Kb=Kp;
Li=0.63E-3;
Ri=0.83;
J_m=4.2E-6;
%% Parametros Helice
R=0.1; % Radio helices
J_T=J_m+1/2*M*R^2;
b=2; % nº de palas
c=1.5/100; % Cuerda media
sigma=b*c/(pi*R);
a=5.73;
theta_0=deg2rad(15);  % 15 grados en raiz
theta_1=deg2rad(-0);
C_d0=0.01;
%% Coeficientes de resistencia aerodinamicos y Areas equivalentes
%         Basicamentes son inventados pero podrian determinarse mediante
%         ensayos en tunel
C_x=1.2;
C_y=C_x;
A_x=0.015;
A_y=A_x;
C_z=2.2;
A_z=0.16;
% Coeficientes de amortiguacion
C_ax=0.1;
C_ay=C_ax;
C_az=0.01;
%% Param generales
rho=1.225;
g=9.81;
v_i0=sqrt( 1/4*m*g/(2*rho*pi*R^2)  ); % Velocidad inducida a punto fijo


