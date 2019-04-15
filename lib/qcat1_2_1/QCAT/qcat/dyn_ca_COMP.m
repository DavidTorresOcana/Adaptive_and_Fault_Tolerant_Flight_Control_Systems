function OUT = dyn_ca_COMP(v,u_ant,W,us,B,plim,rlim,T,Wv,W1,W2,alg, gammma,imax,k,m,tol)
%% Completed Dynamic Control Allocator, using all algorithms of CA
%
% Author: David Torres
% Based on work of Ola Härkegård
% Date: 05/06/2014
% 
% dyn_ca_COMP - Dynamic control allocation 
%
%  [u,W,iter] = dyn_ca_COMP(v,u(t-T),W_ant,us,B,plim,[rlim,T,Wv,W1,W2,alg, gammma,imax,k,m],options)
%
% Performs dynamic control allocation for the actual virtual
% with the previous u(t-T). For each value of v, the control
% signal u is determined by solving:
%     
%   min ||W1(u(t)-us)||^2+ ||W2(u(t)-u(t-T))||^2
%  u in M
%
% where M is the set of control signals solving
%
%   min   ||Wv(Bu-v)||
%  u in U
%
% where U is the set of feasible control signals with respect to
% position and, possibly, rate limits.
%
%  Inputs:
%  -------
% v     commanded virtual control trajectory (k x 1)
% u_ant  Previous optimal control
% us    desired control (m x 1)
% W_ant  Previous active constraints in u
% B     control effectiveness matrix (k x m)
% plim  position limits [min max] (m x 2)
% rlim  rate limits [min max] (m x 2) ([] --> no rate limits)
% T     sampling time [1]
% Wv    virtual control weighting matrix (k x k) [I]
% W1    control position weighting matrix (m x m) [I]
% W2    control rate weighting matrix (m x m) [0]
% alg   Algorithm to use:[ WLS,WLSC,FXP,SLS,MLS,CGI,IP]
% gamma Relative weight of error [1e6]
% imax  Maximum number of iterations
% k and m are the dimensions of the problem B=(k x m)
% tol Tolerance for IP problem
%
%
%  --------
%
%  Outputs:
%  -------
% u     optimal control
%                            0 if u_i not saturated
% Working set syntax: W_i = -1 if u_i = umin_i
%                           +1 if u_i = umax_i
% iter  no. of iterations
%
% See also: QCATDOC ,  QP_SIM , WLS_ALLOC  ,  WLSC_ALLOC  ,  FXP_ALLOC.

if size(B,1)~=k && size(B,2)~=m
    % Retrieve normal Dimensions
    B = reshape(B,k,m);
    plim = reshape(plim,m,2);
    rlim = reshape(rlim,m,2);
    Wv = reshape(Wv,k,k);
    W1 = reshape(W1,m,m);
    W2 = reshape(W2,m,m);
end
  

  % Formulate as a regular QP problem: min ||Wu(u-ud)||
  W1sq = W1^2;
  W2sq = W2^2;
  Wu = sqrtm(W1sq+W2sq);
  invWusq = inv(W1sq+W2sq);

  ud = invWusq*(W1sq*us+W2sq*u_ant);

  % Overall position limits
  if isempty(rlim)
    umin = plim(:,1);
    umax = plim(:,2);
  else
    umin = max(plim(:,1),u_ant+rlim(:,1)*T);
    umax = min(plim(:,2),u_ant+rlim(:,2)*T);
  end
  
  % Use WLS -- fast and robust.
switch lower(alg)
    case 'wls'
        [u,W,iter] = wls_alloc(B,v,umin,umax,Wv,Wu,ud,gammma,u_ant,W,imax);
    case 'wlsc'
        [u,W,iter] = wlsc_alloc(B,v,umin,umax,Wv,Wu,ud,gammma,u_ant,W,imax);
    case 'fxp'
        [u,iter] = fxp_alloc(B,v,umin,umax,Wv,Wu,ud,gammma,u_ant,imax);
        W = 0*ones(size(u));
    case 'sls'
        [u,W,iter] = sls_alloc(B,v,umin,umax,Wv,Wu,ud,u_ant,W,imax); % Gamma should be included in Wv
    case 'mls'
        [u,W,iter] = mls_alloc(B,v,umin,umax,Wv,Wu,ud,u_ant,W,imax); % Gamma should be included in Wv
    case 'cgi'
        [u,iter] = cgi_alloc(B,v,umin,umax,Wv,Wu,ud,imax); % Gamma should be included in Wv
        W = 0*ones(size(u));
    case 'ip'
        [u,iter] = ip_alloc(B,v,umin,umax,ud,gammma,tol,imax);
        W = 0*ones(size(u));
end
    
  OUT = [u ; W ; iter];
  
  
  