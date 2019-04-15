function OUT = qp_ca_COMP(v,u_ant,W,ud,B,plim,rlim,T,Wv,Wu,alg, gammma,imax,k,m,tol)
%% Completed QP Control Allocator, using all algorithms of CA
%
% Author: David Torres
% Based on work of Ola Härkegård
% Date: 05/06/2014
% 
% qp_ca_COMP - QP control allocation 
%
%  [u,W,iter] = qp_ca_COMP(v,u(t-T),W_ant,us,B,plim,[rlim,T,Wv,W1,W2,alg, gammma,imax,k,m],options)
%
% Performs QP control allocation for the actual virtual input
% taking into account rate limits. For each value of v, the control
% signal u is determined by any of the algorithms selected: 
%    [WLS,WLSC,FXP,SLS,MLS,CGI,IP]
% 
% 
%  Inputs:
%  -------
% v     commanded virtual control trajectory (k x 1)
% u_ant  Previous optimal control
% ud    desired control (m x 1)
% W_ant  Previous active constraints in u
% B     control effectiveness matrix (k x m)
% plim  position limits [min max] (m x 2)
% rlim  rate limits [min max] (m x 2) ([] --> no rate limits)
% T     sampling time [1]
% Wv    virtual control weighting matrix (k x k) [I]
% Wu    control weighting matrix (m x m) [I]
% alg   Algorithm to use:  [WLS,WLSC,FXP,SLS,MLS,CGI,IP]
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
% See also: QCATDOC ,  QP_SIM , WLS_ALLOC, WLSC_ALLOC, FXP_ALLOC, SLS_ALLOC, MLS_ALLOC, CGI_ALLOC, IP_ALLOC .
% 

if size(B,1)~=k && size(B,2)~=m
    % Retrieve normal Dimensions
    B = reshape(B,k,m);
    plim = reshape(plim,m,2);
    rlim = reshape(rlim,m,2);
    Wv = reshape(Wv,k,k);
    Wu = reshape(Wu,m,m);

end

  % Formulate as a regular QP problem: min ||Wu(u-ud)||
  
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
  
  
  
  