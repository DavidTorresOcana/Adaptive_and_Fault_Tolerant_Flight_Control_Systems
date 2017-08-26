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

x_k_1= transpose( [0.1,0,0,0,0,0]);  % Este valor puede estar dado o por la integracion o por los sensores
x_k= transpose( [0.1,0,0,0,0,0]);
r_k_1=transpose( [0.5,0,0,0,0,0]);  % Referencia dada
At=0.1;

initial_K=K; % Cojemos por ahora la K del fabricante; La inicial deberia ser la del paso anterior

initial_K = reshape(initial_K,size(initial_K,1)*size(initial_K,2),1)  % Para usar fmingc hay que usar un ector como parametro a optimizar

[J,Grad] = FuncionCoste(A,B,At,x_k_1,x_k,r_k_1,initial_K);

Grad=reshape(Grad,size(K));
display(J)
display(Grad)
fprintf(' Coste inicial y gradiente\n');
fprintf(' Se procedera a encontrar la K optima\n');
% pause
FunciondeCoste = @(t) FuncionCoste(A,B,At,x_k_1,x_k,r_k_1,t);

options = optimset('MaxIter', 1);
[K_optima, Coste] = fmincg(FunciondeCoste, initial_K, options);

fprintf(' K optima y evolucion del coste\n');

K_optima=reshape(K_optima,size(K));
display(K_optima)
display(K)


fprintf(' Pasaremos a hacer una simulacion\n');
pause
%% Simulacion
% Intentaremos integrar con un RK3 las ecuaciones usando ambas K's
x(:,1)= transpose( [0.1,0,0,0,0,0]);  % Posicion inicial
x(:,2)= transpose( [0.1,0,0,0,0,0]);
r= transpose( [0.5,0,0,0,0,0]);

% K del fabricante

for j=1:40
%     r=j/100 +0.1;
    k1=A*x(:,j) - B*K*(x(:,j)-r);
    k2=A*( x(:,j) + At*k1/2)  - B*K*(  x(:,j) + At*k1/2 -r);
    k3=A*( x(:,j) + At*(2*k1-k1) )  - B*K*(  x(:,j) + At*(2*k1-k1) -r);
    x(:,j+1) = x(:,j) +At/6*( k1+4*k2 + k3 );
    
    
end

subplot(2,1,1)
plot(x(1,:))
hold on

% plot(x(4,:),'g')
hold off

fprintf('K discreta optima\n');
% pause

% K Opt discreta

for j=2:40
%     r=j/100 + 0.1;
    
    initial_K=K; % La primera K la cogemos del LQR normal. Hacemos la optimizacion desde la K anterior
    % Se podria evitar oscilaciones de la solucion optimizando siempre
    % desde una misma K_opt por ejemplo la K del LQR normal
    
    initial_K = reshape(initial_K,size(initial_K,1)*size(initial_K,2),1);

% % %     FunciondeCoste = @(t) FuncionCoste(A,B,At,x(:,j-1),x(:,j),r,t);
% % %     [K_optima, Coste] = fmincg(FunciondeCoste, initial_K, options);

    K_optima = LQRDiscretoFUNC(At,x(:,j-1),x(:,j),r,initial_K);
    
    
    K=reshape(K_optima,size(K));
    
    
    k1=A*x(:,j) - B*K*(x(:,j)-r);
    k2=A*( x(:,j) + At*k1/2)  - B*K*(  x(:,j) + At*k1/2 -r);
    k3=A*( x(:,j) + At*(2*k1-k1) )  - B*K*(  x(:,j) + At*(2*k1-k1) -r);
    x(:,j+1) = x(:,j) +At/6*( k1+4*k2 + k3 );

end

subplot(2,1,2)
plot(x(1,:))
hold on

% plot(x(4,:),'g')
hold off

end

function [J,Grad]=FuncionCoste(A,B,At,x_k_1,x_k,r_k_1,K)

    K=reshape(K, 4 , 6); % CUIDADO CON ESTO!!!!  Cuando se cambien las dimesiones de los vectores
    
    % Se observa que puede tener un mejor comportamiento cuando el estado
    % actual (K) se desprecia y se estima a partir del antrior x_k_1
% % % % % % %     x_k=( eye(6)+At*A )*x_k_1 -At*B*K*(  x_k_1 - r_k_1  );  %
% Cambiar esto cuando se queira!!!! 
    
    J= transpose (  ( eye(6)+At*A )*x_k -At*B*K*(  x_k - r_k_1  )  -r_k_1 ) * (  ( eye(6)+At*A )*x_k -At*B*K*(  x_k - r_k_1  ) -r_k_1 );
    
    
% %     Grad=-2*(  ( eye(6)+At*A )*x_k_1 -At*B*K*(x_k_1-r_k_1) - r_k_1   )*At*B* [ (x_k_1-r_k_1)'; (x_k_1-r_k_1)' ; (x_k_1-r_k_1)' ; (x_k_1-r_k_1)' ]  ;
    for i=1:size(B,2)
        
        Grad ( i, : ) = 2*(  ( eye(6)+At*A )*x_k - At*B*K*(  x_k - r_k_1  ) - r_k_1 )' * (  -At*B(:,i)*( x_k -r_k_1 )'  - ( ( eye(6)+At*A )-At*B*K ) * At*B(:,i)* (x_k_1 -r_k_1 )'  );
        
    end
    
    Grad = reshape(Grad,size(Grad,1)*size(Grad,2),1);  % Para usar fmingc hay que usar un vector como gradiente tambien
    
end