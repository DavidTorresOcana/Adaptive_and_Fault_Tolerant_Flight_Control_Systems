function LQRDiscretoControlador()
clc;
clear all
global A B
% Set the model parameters of the 3DOF HOVER.
% These parameters are used for model representation and controller design.
[ Ktn, Ktc, Kf, l, Jy, Jp, Jr, g ] = setup_hover_configuration();
%
% For the following state vector: X = [ theta; psi; theta_dot; psi_dot]
% Initialization the state-Space representation of the open-loop System
HOVER_ABCD_eqns;
%% LQR suministrado por el fabricante
    Q = diag([500 350 350 0 20 20] );
    R = 0.01*diag([1 1 1 1]);
    % Automatically calculate the LQR controller gain
    K = lqr( A, B, Q, R )    
    
    
    
%% LQR calculado optimizando coste de forma discreta
% Definimos el vector de estado x =[psi,theta,phi,psi',theta',phi']

x_k= zeros(6,1);  % Este valor puede estar dado o por la integracion o por los sensores
r_k=transpose( [0,0,0,0.1,0,0]);  % Referencia dada en el paso anterior
r_k1= transpose( [0,0,0,0.5,0,0]); % Referencia a donde queremos ir
At=0.01;

initial_K=K; % Cojemos por ahora la K del fabricante; La inicial deberia ser la del paso anterior

initial_K = reshape(initial_K,size(initial_K,1)*size(initial_K,2),1)  % Para usar fmingc hay que usar un ector como parametro a optimizar

[J,Grad] = FuncionCoste(A,B,At,x_k,r_k,r_k1,initial_K);

Grad=reshape(Grad,size(K));
display(J)
display(Grad)
fprintf(' Coste inicial y gradiente\n');
fprintf(' Se procedera a encontrar la K optima\n');
% pause
FunciondeCoste = @(t) FuncionCoste(A,B,At,x_k,r_k,r_k1,t);

options = optimset('MaxIter', 100);
[K_optima, Coste] = fmincg(FunciondeCoste, initial_K, options);

fprintf(' K optima y evolucion del coste\n');

K_optima=reshape(K_optima,size(K));
display(K_optima)
display(K)
plot(Coste)

close all
fprintf(' Pasaremos a hacer una simulacion\n');
pause
%% Simulacion
% Intentaremos integrar con un RK3 las ecuaciones usando ambas K's
x(:,1)= transpose( [0,0,0,0,0.1,0]);  % Posicion inicial
r=transpose( [0,0,0,0,0.1,0]);  % Referencia dada en el primer paso
r_obj= transpose( [0,0,0,0,0.5,0]);

% K del fabricante
for j=1:20
    if j>1 
        r=r_obj;
    end

    k1=A*x(:,j) - B*K*(x(:,j)-r);
    k2=A*( x(:,j) + At*k1/2)  - B*K*(  x(:,j) + At*k1/2 -r);
    k3=A*( x(:,j) + At*(2*k1-k1) )  - B*K*(  x(:,j) + At*(2*k1-k1) -r);
    x(:,j+1) = x(:,j) +At/6*( k1+4*k2 + k3 );
    
    
end
subplot(2,1,1)
plot(x(2,:))
hold on

plot(x(5,:),'g')
hold off

fprintf('K discreta optima\n');
% pause

% K Opt discreta

for j=1:20
    if j>1 
        r=r_obj;
    end
    initial_K=K; % La primera K la cogemos del LQR normal. Hacemos la optimizacion desde la K anterior
    % Se podria evitar oscilaciones de la solucion optimizando siempre
    % desde una misma K_opt por ejemplo la K del LQR normal
    
    initial_K = reshape(initial_K,size(initial_K,1)*size(initial_K,2),1);
  
% % % % % %     FunciondeCoste = @(t) FuncionCoste(A,B,At,x(:,j),r,r_obj,t);
% % % % % %     [K_optima, Coste] = fmincg(FunciondeCoste, initial_K, options);

    K_optima = LQRDiscretoFUNC(At,x(:,j),r,r_obj,initial_K);
    
    
    K=reshape(K_optima,size(K));
    
    
    k1=A*x(:,j) - B*K*(x(:,j)-r);
    k2=A*( x(:,j) + At*k1/2)  - B*K*(  x(:,j) + At*k1/2 -r);
    k3=A*( x(:,j) + At*(2*k1-k1) )  - B*K*(  x(:,j) + At*(2*k1-k1) -r);
    x(:,j+1) = x(:,j) +At/6*( k1+4*k2 + k3 );

end
subplot(2,1,2)
plot(x(2,:))
hold on

plot(x(5,:),'g')
hold off

end
function [J,Grad]=FuncionCoste(A,B,At,x_k,r_k,r_k1,K)

    K=reshape(K, 4 , 6); % CUIDADO CON ESTO!!!!  Cuando se cambien las dimesiones de los vectores
    J= transpose (  ( eye(6)+At*A )*x_k -At*B*K*(x_k-r_k) - r_k1  )* (  ( eye(6)+At*A )*x_k -At*B*K*(x_k-r_k) - r_k1    );
    
    
%     Grad=-2*(  ( eye(6)+At*A )*x_k -At*B*K*(x_k-r_k) - r_k1   )*At*B* [ (x_k-r_k)'; (x_k-r_k)' ; (x_k-r_k)' ; (x_k-r_k)' ]  ;
    for i=1:size(B,2)
        Grad ( i, : ) = -2*(  ( eye(6)+At*A )*x_k -At*B*K*(x_k-r_k) - r_k1   )' * At*B(:,i) * (x_k-r_k)';
    end
    
    Grad = reshape(Grad,size(Grad,1)*size(Grad,2),1);  % Para usar fmingc hay que usar un vector como gradiente tambien
    
end