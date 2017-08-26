function [ VEFI ] = VE_to_VEFI( VE)

p=VE(1);
q=VE(2);
r=VE(3);
q1=VE(4);
q2=VE(5);
q3=VE(6);
q4=VE(7);

%% Actitud en Angulos de Euler y sus derivadas. 
qp = 1/2*HamiltonProduct( [q1;q2;q3;q4],[0;p;q;r] ); 

VEFI(1) = atan2(2*q3*q4+2*q1*q2,2*q1^2+2*q4^2-1 ); % Balance phi
VEFI(2) = -asin( 2*q2*q4-2*q1*q3 ); % Angulo asiento theta
VEFI(3) = atan2( 2*q2*q3+2*q1*q4,2*q1^2+2*q2^2-1 ); % Angulo de guiñada psi


VEFI(4) = 2/(1+(   (2*q3*q4+2*q1*q2)/(2*q1^2+2*q4^2-1)  )^2)*[ q2*(q1^2-2*q4^2+1)+2*q1*q3*q4  ,  q1*(2*q1^2+2*q4^2-1)  ,  q4*(2*q1^2+2*q4^2-1)  ,  -(-q1^2*q3 +4*q1*q2*q4+2*q3*q4^2+q3 ) ]*qp* 1/(2*q1^2+2*q4^2-1)^2; % Derivada de phi
VEFI(5) = -2/( sqrt(1-4*(q2*q4-q1*q3)^2)  )*[-q3,q4,-q1,q2]*qp;  % Theta prima
VEFI(6) = 2/(1+(  (2*q2*q3+2*q1*q4)/(2*q1^2+2*q2^2-1)  )^2)*[  -(2*q1^2*q4+4*q1*q2*q3-4*q2*q4+q4  ) , q3*(2*q1^2-1)-4*q1*q4  ,  q2  , q1 ]*qp*  1/(2*q1^2+2*q2^2-1)^2;  % Psi prima


end