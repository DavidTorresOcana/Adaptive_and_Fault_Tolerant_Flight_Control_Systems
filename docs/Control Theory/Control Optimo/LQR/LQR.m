clear all
close all
clc

load('pss6');
SYSLIN=pss6;
[A,B,C,D,K] = idssdata(SYSLIN);
% pause
%% Sistema sin control

if rank(ctrb(SYSLIN))==size(SYSLIN(1,1).a,1)
    fprintf(' El sistema es controlable!!!!')
else 
    fprintf(' El sistema NO es controlable!!!!')
end
% pause
SYSLIN=ss(SYSLIN)
tf(SYSLIN)
% pause
step(SYSLIN)
% pause
%% Diseño LQR
KMIO=lqr(SYSLIN,diag([1000,1000,200,0,0,0]),0.01*eye(4,4))
pause

% Comparacion
fprintf('  Polos sistema sin control');
eig(SYSLIN)
pause
fprintf('  Polos sistema con LQR');
eig(A-B*KMIO)
pause

% Plots
subplot(1,2,1)
plot(eig(SYSLIN),'*')
title(' Polos sistema sin control')
subplot(1,2,2)
plot(eig(A-B*KMIO),'*')
title(' Polos sistema con LQR')