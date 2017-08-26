function LQRDiscretoGradiente()
clc;
clear all
%%Comentarios de este metodo
% Se confia solo un parametro la accion de ponderacion entre la K calculada
% y el Gradiente obtenido de la prediccion d ela accion.
% El dejar solo un aprmetro de ponderacion seria MENOS efectivo usar los
% parametro Q y R en el LQR que son parametro de ponderacion. Lo ideal
% seria poder calcularlos en funcion de X_actual y r_requerida
[ Ktn, Ktc, Kf, l, Jy, Jp, Jr, g ] = setup_hover_configuration();

HOVER_ABCD_eqns;
%% LQR suministrado por el fabricante
    Q = diag([500 350 350 0 20 20] );
    R = 0.01*diag([1 1 1 1]);
    % Automatically calculate the LQR controller gain
    K = lqr( A, B, Q, R )    
    
    
    
%% LQR calculado optimizando coste de forma discreta
% Definimos el vector de estado x =[psi,theta,phi,psi',theta',phi']

x_k= zeros(6,1);  % Este valor puede estar dado o por la integracion o por los sensores
r_k=transpose( [0,0,0,0,0,0.1]);  % Referencia dada en el paso anterior
r_k1= transpose( [0,0,0,0,0,0.5]); % Referencia a donde queremos ir
At=0.01;
beta=10;% Se confia a beta la accion de ponderar la importancia del gradiente


initial_K=K; % Cojemos por ahora la K del fabricante; La inicial deberia ser la del paso anterior

initial_K = reshape(initial_K,size(initial_K,1)*size(initial_K,2),1);  

[J,Grad] = FuncionCoste(A,B,At,x_k,r_k,r_k1,initial_K);

Grad=reshape(Grad,size(K));
display(J)
display(Grad)
fprintf(' Coste inicial y gradiente\n');
fprintf(' Se procedera a encontrar la K optima\n');
% pause

K_optima=reshape(initial_K,size(K)) + beta*reshape(Grad,size(K)); % Se confia a beta la accion de ponderar la importancia del gradiente

fprintf(' K optima y evolucion del coste\n');

display(K_optima)
display(K)

close all
fprintf(' Pasaremos a hacer una simulacion\n');
pause
%% Simulacion
% Intentaremos integrar con un RK3 las ecuaciones usando ambas K's
x(:,1)= zeros(6,1);  % Psicion inicial
r=transpose( [0,0,0.1,0,0,0.1]);  % Referencia dada en el primer paso
r_obj= transpose( [0,0,0.5,0,0,0.5]);

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
plot(x(6,:))
hold on
plot(x(2,:),'r')
plot(x(3,:),'g')
hold off

fprintf('K discreta optima\n');
% pause

% K Opt discreta
for j=1:20
    if j>1 
        r=r_obj;
    end
    
    K = reshape(K,size(K,1)*size(K,2),1);
    [Coste,Grad] = FuncionCoste(A,B,At,x(:,j),r,r_obj,K);
    
    K =  reshape(K,4,6) + beta*reshape(Grad,4,6)
    
    k1=A*x(:,j) - B*K*(x(:,j)-r);
    k2=A*( x(:,j) + At*k1/2)  - B*K*(  x(:,j) + At*k1/2 -r);
    k3=A*( x(:,j) + At*(2*k1-k1) )  - B*K*(  x(:,j) + At*(2*k1-k1) -r);
    x(:,j+1) = x(:,j) +At/6*( k1+4*k2 + k3 );

end
subplot(2,1,2)
plot(x(6,:))
hold on
plot(x(2,:),'r')
plot(x(3,:),'g')
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