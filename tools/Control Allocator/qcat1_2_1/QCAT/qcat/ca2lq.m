function R1 = ca2lq(B,R2,Wu)
  
% CA2LQ - Incorporate control allocator into LQ controller.
% 
%  R1 = ca2lq(B,R2,Wu)
% 
% Given a dynamic system
%  .                             .
%  x = Ax + B1u, B1 = B2*B  <=>  x = Ax + B2v, v = Bu
%
% this function calculates the weighting matrix R1 such that the
% design criterium
% 
%  min Integral (x'Qx + u'R1u) dt   (LQ control)
%   u
%
% gives the same linear optimal control law as
% 
%  min Integral (x'Qx + v'R2v) dt   (LQ control)
%   v
% 
%  min ||Wu u||  subj. to  Bu = v   (control allocation)
%   u
%
% Controller structures:
%
%    x   -------  u           x   -------  v   ----  u
%   --->| Q, R1 |--->   <=>  --->| Q, R2 |--->| Wu |--->
%        -------                  -------      ----
%        LQ ctl                   LQ ctl     Ctl alloc
%
% See also LQ2CA, LQR, WPINV.
  
% Thesis, Theorem 10.4:
  R1 = Wu^2+B'*(R2-inv(B*inv(Wu)^2*B'))*B;