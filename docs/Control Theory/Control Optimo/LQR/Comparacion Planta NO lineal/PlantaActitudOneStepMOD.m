function [ SALIDA ] = PlantaActitudOneStepMOD( VE_ant,VEFI_ant,U,omega,AT )
% PlantaCompletaOneStep es la funcion que ajecuta un paso de integracion
% con un RK4 del modelo de quadcopter completo desde un VE_ant en t hasta
% un t + AT  obteniendo VE
% Los parametros del sistema se puden cambiar y ver en la funcion
% SetupCompleto()
%   VE= Vector estado completo En ejes cuerpo incluyendo cuaternion. 
%   U = Señal de control (V)
%         

%         | p |
%         | q |
%   VE =  | r |_
%         | q1  |
%         | q2  |
%         | q3  |_
%         | q4 |

%         |
%         | I^(-1)*(  -1/2*M*R^2*[0;0; omegap_2+omegap_4-omegap_1-omegap_3] -[q*r*(Iz-Iy);pr*(Ix-Iz);pq*(Iy-Ix)] - 1/2*M*R^2*[q;-p;0]*(omega_2+omega_4-omega_1-omega_3)+ M_g+ (M_Cf_a_1+M_Cf_a_2+M_Cf_a_3+M_Cf_a_4)+(M_Cm_a_1+M_Cm_a_2+M_Cm_a_3+M_Cm_a_4)   )
%d VE/dt =|                                                                                                                                                     
%         |                                                 |
%         | 1/2*HamiltonProduct( [q1;q2;q3;q4],[0;p;q;r] )   |
%         |                                                   |
%         |                                                  |

global  v_i0 options m M Ix Iy Iz l h Kp Kb Li Ri R b c sigma a theta_0 theta_1 C_d0 A_x A_y A_z C_x C_y C_z g rho  J_T r_cg C_ax C_ay C_az

%% Evaluamos los coeficientes aerodinamicos (obteniendo fuerzas y momentos)
% e integramos la dinamica del motor/helice  
r1=[l;0;-h];
r2=[0;l;-h];
r3=[-l;0;-h];
r4=[0;-l;-h];

[omegaNEXT(1),F_a(:,1),M_Cf_a(:,1),M_Cm_a(:,1),omegap(1) ] =Rotor_i(r1,U(1),omega(1),VE_ant,VEFI_ant(3),AT);
[omegaNEXT(2),F_a(:,2),M_Cf_a(:,2),M_Cm_a(:,2),omegap(2) ] =Rotor_i(r2,U(2),omega(2),VE_ant,VEFI_ant(3),AT);
[omegaNEXT(3),F_a(:,3),M_Cf_a(:,3),M_Cm_a(:,3),omegap(3) ] =Rotor_i(r3,U(3),omega(3),VE_ant,VEFI_ant(3),AT);
[omegaNEXT(4),F_a(:,4),M_Cf_a(:,4),M_Cm_a(:,4),omegap(4) ] =Rotor_i(r4,U(4),omega(4),VE_ant,VEFI_ant(3),AT);
% disp(M_Cm_a)
% F_a
% M_Cf_a
% M_Cm_a
%% Ejecutamos un paso de integracion con RK4 
k1= F ( VE_ant ,omegap,F_a,M_Cf_a,M_Cm_a,omega );

k2= F ( VE_ant +1/2*AT*k1' ,omegap,F_a,M_Cf_a,M_Cm_a,omega );
k3= F ( VE_ant +1/2*AT*k2' ,omegap,F_a,M_Cf_a,M_Cm_a,omega );
k4= F ( VE_ant +AT*k3' ,omegap,F_a,M_Cf_a,M_Cm_a,omega );

VE  = VE_ant + AT/6* ( k1' + 2*k2' + 2*k3' +k4');
% Normalizacion del cuaternion
VE(4:7)=VE(4:7)./norm(VE(4:7));

%% Evaluamos el vector estado en Ejes tierra (F_I)
% Una vez se ha obtenido el Nuevo VE, se obtine el VEFI en ejes tierra
% (Inerciales)

VEFI= VE_to_VEFI( VE );
SALIDA=[VE;VEFI';omegaNEXT'];


end


function [ DIFFVE ] = F( VE_ant,omegap,F_a,M_Cf_a,M_Cm_a,omega )
% Hacemos una funcion que evalua la derivada del vector estado ( d VE/dt )
% ya que hay que evaluar el mismo en diferentes valores devido al RK4. Es
% la forma mas comoda
global  v_i0 options m M Ix Iy Iz l h Kp Kb Li Ri R b c sigma a theta_0 theta_1 C_d0 A_x A_y A_z C_x C_y C_z g rho  J_T r_cg C_ax C_ay C_az AT

p=VE_ant(1);
q=VE_ant(2);
r=VE_ant(3);
q1=VE_ant(4);
q2=VE_ant(5);
q3=VE_ant(6);
q4=VE_ant(7);

F_g=HamiltonProduct([q1;-q2;-q3;-q4],HamiltonProduct([0;0;0;g],[q1;q2;q3;q4]) ); % q* . g .q
M_g =cross( r_cg, F_g(2:4) );
DIFFVE(1:3) =[1/Ix,0,0;0,1/Iy,0;0,0,1/Iz]*(  -1/2*M*R^2*[0;0; omegap(2)+omegap(4)-omegap(1)-omegap(3)]- [C_ax*p;C_ay*q;C_az*r] -[q*r*(Iz-Iy);p*r*(Ix-Iz);p*q*(Iy-Ix)] - 1/2*M*R^2*[q;-p;0]*(omega(2)+omega(4)-omega(1)-omega(3))+ M_g + (M_Cf_a(:,1)+M_Cf_a(:,2)+M_Cf_a(:,3)+M_Cf_a(:,4))+   (M_Cm_a(:,1)+M_Cm_a(:,2)+M_Cm_a(:,3)+M_Cm_a(:,4))   );
DIFFVE(4:7) =1/2*HamiltonProduct( [q1;q2;q3;q4],[0;p;q;r] );
end

function [omegaNEXT,F_a,M_Cf_a,M_Cm_a,omegap ] =Rotor_i(ri,U,omega,VE_ant,z_FI,AT)
% Daría igual integrar la dinamica del rotor con un esquema Euler pero
% somos burros y utilizamos un RK4. Su dinamica es muy rapida!!!
global  v_i0 options m M Ix Iy Iz l h Kp Kb Li Ri R b c sigma a theta_0 theta_1 C_d0 A_x A_y A_z C_x C_y C_z g rho  J_T r_cg C_ax C_ay C_az

p=VE_ant(1);
q=VE_ant(2);
r=VE_ant(3);

%Velocidades
V_A=  cross([p;q;r],ri);
mu = norm([V_A(1);V_A(2)])/(omega*R);
v_i = ModeloAB( norm([V_A(1);V_A(2)]) , - V_A(3) ); % Hay que meter la velocidad vertical como positiva hacia arriba??
lambda_i =  v_i/(omega*R);
lambda = lambda_i - V_A(3)/(omega*R);

% Coeficientes
C_T=sigma*a/4*(theta_0*(2/3+mu^2) +theta_1/2*(1+mu^2) +  lambda );
% % mux=mu; muy=0; ???
% % C_T=-2*0.745*(-0.15)*sqrt( mux^2+muy^2+0.447^2*( V_A(3)/(omega*R) )^2 + lambda^2 )
% % C_T=2*v_i/(omega*R)*sqrt( mu^2+(v_i/(omega*R) - V_A(3)/(omega*R))^2  )
C_Q=    +sigma*a/4*( 2/3*theta_0 + lambda + theta_1/2)*lambda + sigma*C_d0/8*(1+mu^2);   % Positivo siempre!!! 
% % C_Q= -lambda*C_T
% % pause
C_H= -sigma*a/4*mu*lambda*(theta_0+theta_1/2) + sigma*C_d0/4*mu;
C_Mh= sigma*a/4*( 2/3*theta_0*mu+1/2*lambda*mu+theta_1*mu/2  );
%% Dada un V , w_h y omega calculamos los coeficientes de la helice
% Fuerzas
%                                |-C_H*[V_A(1);V_A(2)]/norm([V_A(1);V_A(2)])  |
% F_a =   rho*pi*R^2*(omega*R)^2*|                                            | 
%                                |-1/(1-R^2/(16*z_FI^2))*C_T                  |  
if V_A(1)==0 && V_A(2)==0
    H=[0;0];
else
    H=C_H*(-[V_A(1);V_A(2)])./norm([V_A(1);V_A(2)]);
end

F_a =   rho*pi*R^2*(omega*R)^2*   [  H ; -C_T ];
% pause
% Momentos
if ri(2)==0   % Discernimos el sentido de giro de la helice
   signo=1;
else
    signo=-1;
end
if V_A(1)==0 && V_A(2)==0
    Mh=[0;0];
else
    Mh=signo*C_Mh*(-[V_A(1);V_A(2)])./norm([V_A(1);V_A(2)]);% Nota!!!! Depende el signo de C_Mh del sentido de giro???
end
M_Cm_a = rho*pi*R^3*(omega*R)^2* [  Mh  ;  signo*C_Q ] ;

M_Cf_a = cross( ri, F_a ); % Momentos de fuerzas!!!!!
                     
%% Dinamica del motor/helice 
        % Modelo Semicompleto
        
        omegap= 1/(Ri/Kp*J_T+Li/Kp*2*rho*pi*R^2*C_Q*(omega))*(U-Ri/Kp*rho*pi*R^5*C_Q*(omega)^2-Kb*(omega)  );
        %  RK4    
        k1= omegap;
        
        k2 = 1/(Ri/Kp*J_T+Li/Kp*2*rho*pi*R^2*C_Q*(omega+1/2*k1*AT))*(U-Ri/Kp*rho*pi*R^5*C_Q*(omega+1/2*k1*AT)^2-Kb*(omega+1/2*k1*AT)  );
        k3 = 1/(Ri/Kp*J_T+Li/Kp*2*rho*pi*R^2*C_Q*(omega+1/2*k2*AT))*(U-Ri/Kp*rho*pi*R^5*C_Q*(omega+1/2*k2*AT)^2-Kb*(omega+1/2*k2*AT)  );
        k4 = 1/(Ri/Kp*J_T+Li/Kp*2*rho*pi*R^2*C_Q*(omega+k3*AT))*(U-Ri/Kp*rho*pi*R^5*C_Q*(omega+k3*AT)^2-Kb*(omega+k3*AT)  );
        % Modelo Simplificado (Quizá demasiado simple y rapido)
%         omegap= Kp/(2*Li)*( (U/omega-Kb)/(rho*pi*R^5*C_Q) - Ri/Kp*omega   );
%         %  RK4    
%         k1= omegap;
%         k2= Kp/(2*Li)*( (U/(omega+1/2*AT*k1)-Kb)/(rho*pi*R^5*C_Q) - Ri/Kp*(omega+1/2*AT*k1)  );
%         k3= Kp/(2*Li)*( (U/(omega+1/2*AT*k2)-Kb)/(rho*pi*R^5*C_Q) - Ri/Kp*(omega+1/2*AT*k2)  );
%         k4= Kp/(2*Li)*( (U/(omega+AT*k3)-Kb)/(rho*pi*R^5*C_Q) - Ri/Kp*(omega+AT*k3)  );
omegaNEXT  = omega + AT/6* ( k1 + 2*k2 + 2*k3 +k4);
% %  Euler
%  omegaNEXT  = omega + AT* omegap;
% omegap
% omega
% omegaNEXT
% pause
if C_T<0 || C_Q<0 || lambda<0 
    fprintf('eeeeinnng');
    pause;
end

end
function v_i=ModeloAB( Vx , Vz )
% Resolucion numerica del modelo de Velocidad inducidad AB
global v_i0 options

    v_i = fsolve(@FuncionAB,v_i0,options);
%     v_i = fsolve(@FuncionVI,v_i0,options);
    
    function FF=FuncionAB( v_i )
        FF = 0.745*v_i/v_i0* sqrt( (Vx/v_i0)^2 + 0.447^2*(Vz/v_i0)^2+(Vz/v_i0+v_i/v_i0 )^2    ) -1;
    end
%     function FF=FuncionVI( v_i )
%         FF = v_i/v_i0* sqrt( (Vx/v_i0)^2 + (Vz/v_i0)^2+(Vz/v_i0+v_i/v_i0 )^2    ) -1;
%     end


% 
% Vi2=v_i/v_i0
% Vz2=Vz/v_i0
% 
% % pause(0.2)

end
