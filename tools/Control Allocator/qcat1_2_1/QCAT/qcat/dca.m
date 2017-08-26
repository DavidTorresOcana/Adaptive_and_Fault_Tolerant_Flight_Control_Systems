function Gvu = dca(B,S,W1,W2,T)
  
% DCA - Compute dynamic control allocation filter.
%
%  Gvu = dca(B,S,W1,W2,[T])
% 
% Computes the explicit filter solution u(t) = Gvu(q)v(t) to the
% unconstrained dynamic control allocation problem
%
%   min  ||W1(u(t)-Sv(t))||^2 + ||W2(u(t)-u(t-T))||^2
%   u(t)
%
%   subj. to  Bu(t)=v(t)
%
%  Inputs:
%  -------
% B     control effectiveness matrix (k x m)
% S     desired steady state control distribution (m x k)
% W1    control position weighting matrix (m x m)
% W2    control rate weighting matrix (m x m)
% T     sampling time [-1]
%
%  Output:
%  ------
% Gvu   optimal filter: u(t) = Gvu(q)*v(t)
% 
%  Example:
%
%   Gvu = dca([1 1 1],[1 0 0]',diag([0 1 1]),diag([3 2 0]),.1)
%   figure(1),step(Gvu)
%   figure(2),bodemag(Gvu)  
%
% See also: DYN_SIM.
  
  if nargin<5,
    T = -1;
  end

  % Number of control inputs.
  m = size(B,2);
  % Net weighting matrix.
  W = sqrtm(W1^2+W2^2);
  Wi = inv(W);
  % u(t) = ESv(t)+Fu(t-T)+Gv(t)
  G = Wi*pinv(B*Wi);
  E = (eye(m)-G*B)*Wi^2*W1^2;
  F = (eye(m)-G*B)*Wi^2*W2^2;
  % Control allocation transfer function.
  Gvu = ss(F,F*(E*S+G),eye(m),E*S+G,T);
  Gvu = tf(Gvu);
