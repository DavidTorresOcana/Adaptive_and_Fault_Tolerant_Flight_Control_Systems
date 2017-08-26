function [q]=QuaternionfromEuler(phi,theta,psi)
% options=optimset('MaxFunEvals',3000);
% q=fsolve(@QUATERNIONEULER,[0.6,0.4,0.5,0.4]); % Hay que acertar con los angulos poniendo diferentes valores iniciales del quaternion
%     function F=QUATERNIONEULER(q)
%         F=[ tan(phi)-(2*q(3)*q(4)-2*q(1)*q(2))/(2*q(1)^2+2*q(4)^2-1)  ;
%             sin(-theta)-2*q(2)*q(4)-2*q(1)*q(3);
%             tan(psi)-(2*q(2)*q(3)-2*q(1)*q(4))/(2*q(1)^2+2*q(2)^2-1);
%             q(1)^2+q(2)^2+q(3)^2+q(4)^2-1];
%     end  % Este algoritmo da q* !!!!!

q=[ cos(psi/2)*cos(theta/2)*cos(phi/2)+ sin(psi/2)*sin(theta/2)*sin(phi/2) ;
      cos(psi/2)*cos(theta/2)*sin(phi/2)- sin(psi/2)*sin(theta/2)*cos(phi/2) ;
       cos(psi/2)*sin(theta/2)*cos(phi/2)+ sin(psi/2)*cos(theta/2)*sin(phi/2)   ;
        sin(psi/2)*cos(theta/2)*cos(phi/2)- cos(psi/2)*sin(theta/2)*sin(phi/2)      ]; % Esto nos da q !!!
end