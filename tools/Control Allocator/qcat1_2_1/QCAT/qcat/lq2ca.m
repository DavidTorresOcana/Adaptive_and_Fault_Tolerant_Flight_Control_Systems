function [R2,Wu] = lq2ca(B,R1)
  
% LQ2CA - Extract control allocator from a given LQ controller.
% 
%  [R2,Wu] = lq2ca(B,R1)
% 
% Given a dynamic system
%  .                             .
%  x = Ax + B1u, B1 = B2*B  <=>  x = Ax + B2v, v = Bu
%
% this function calculates weighting matrices R2 and Wu such that the
% design criterium
% 
%  min Integral (x'Qx + v'R2v) dt   (LQ control)
%   v
% 
%  min ||Wu u||  subj. to  Bu = v   (control allocation)
%   u
% 
% gives the same linear optimal control law as
% 
%  min Integral (x'Qx + u'R1u) dt   (LQ control)
%   u
%
% Controller structures:
%
%    x   -------  u           x   -------  v   ----  u
%   --->| Q, R1 |--->   <=>  --->| Q, R2 |--->| Wu |--->
%        -------                  -------      ----
%        LQ ctl                   LQ ctl     Ctl alloc
%
% See also CA2LQ, LQR, WPINV.
  
% Thesis, Theorem 10.4:
  R2 = inv(B*inv(R1)*B');
  Wu = sqrtm(R1);
  