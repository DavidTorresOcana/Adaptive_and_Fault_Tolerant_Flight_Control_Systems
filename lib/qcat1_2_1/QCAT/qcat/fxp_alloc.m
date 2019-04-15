function [u,iter] = fxp_alloc(B,v,umin,umax,Wv,Wu,ud,gam,u,imax)

% FXP_ALLOC - Control allocation using a fixed-point iterations.
%
%  [u,iter] = fxp_alloc(B,v,umin,umax,[Wv,Wu,ud,gamma,u0,imax])
%
% Solves the weighted, bounded least-squares problem
% 
%   min ||Wu(u-ud)||^2 + gamma ||Wv(Bu-v)||^2
%
%   subj. to  umin <= u <= umax
% 
% using a fixed-point iteration algorithm.
%
%  Inputs:
%  -------
% B     control effectiveness matrix (k x m)
% v     commanded virtual control (k x 1)
% umin  lower position limits (m x 1)
% umax  upper position limits (m x 1)
% Wv    virtual control weighting matrix (k x k) [I]
% Wu    control weighting matrix (m x m) [I]
% ud    desired control (m x 1) [0]
% gamma weight (scalar) [1e3]
% u0    initial point (m x 1)
% imax  no. of iterations [100]
% 
%  Outputs:
%  -------
% u     (approximately) optimal solution
% iter  no. of iterations
%
% See also: WLS_ALLOC, WLSC_ALLOC, IP_ALLOC, QP_SIM.
  
% Number of variables
  m = length(umin);
  
  % Set default values of optional arguments
  if nargin < 10
    imax = 100; % Heuristic value
    [k,m] = size(B);
    if nargin < 9, u = (umin+umax)/2; end
    if nargin < 8, gam = 1e3;         end
    if nargin < 7, ud = zeros(m,1);  end
    if nargin < 6, Wu = eye(m);      end
    if nargin < 5, Wv = eye(k);      end
  end

  % Change of variables.
  x = u - ud;
  xmin = umin - ud;
  xmax = umax - ud;
  
  % Compute the weight epsilon used by Burken et al.
  e = 1/(gam+1);
  
  % m = number of variables (actuators)
  m = size(B,2);
  
  % Notational differences: Q1 = Wu'Wu
  %			    Q2 = Wv'Wv
  %			    CzBu*-dz = v
  
  % Common factor.
  BzT_Q1 = B'*Wv'*Wv;
  % Just above eq. (25).
  H = (1-e) * (BzT_Q1 * B) + e*Wu'*Wu;
  % Eq. (25).
  w = 1/norm(H,'fro');
  % Matrices used for updating solution: 
  % u[k] = sat(Fv-Gu[k-1])
  Fv = (1-e) * w * BzT_Q1 * v;
  G  = w*H - eye(m);
    
  for iter = 1:imax
    % Compute next point without considering the constraints.
    x = Fv - G*x;
    % Saturate to box constraints.
    x = max(xmin,min(xmax,x));
  end
    
  % Back to original variables.
  u = x + ud;
