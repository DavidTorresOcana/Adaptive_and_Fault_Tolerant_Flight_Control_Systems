function [u,iter] = cgi_alloc(B,v,umin,umax,Wv,Wu,ud,imax)

% CGI_ALLOC - Control allocation based on cascading generalized inverses.
%
%  [u,iter] = cgi_alloc(B,v,umin,umax,[Wv,Wu,ud,imax])
%
% Approximately solves the bounded sequential least-squares problem
%
%   min ||Wu(u-ud)||   subj. to   u in M
%
% where M is the set of control signals solving
%
%   min ||Wv(Bu-v)||   subj. to   umin <= u <= umax
%
% using the method of cascading generalized inverses.
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
% imax  max no. of iterations [100]
%
%  Outputs:
%  --------
% u     (approximately) optimal solution
% iter  no. of iterations (= no. of pseudoinverse solutions computed)
%
% See also: SLS_ALLOC, MLS_ALLOC, QP_SIM.
  
% Set default values of optional arguments
  if nargin < 8
    imax = 100; % Heuristic value
    [k,m] = size(B);
    if nargin < 7, ud = zeros(m,1); end
    if nargin < 6, Wu = eye(m);     end
    if nargin < 5, Wv = eye(k);     end
  end

  iter = 1;
  % Compute the optimal solution disregarding the inequality
  % constraints.
  A = Wv*B/Wu;
  b = Wv*(v-B*ud);
  u = ud + Wu\pinv_sol(A,b);

  % Find indeces of infeasible variables.
  i_min = u < umin;
  i_max = u > umax;
  % Set these variables to their limits.
  u(i_min) = umin(i_min);
  u(i_max) = umax(i_max);
  % Let the remaining variables be free.
  i_free = ~(i_min | i_max);
  
  % If the preceeding pseudoinverse solution yielded some variables
  % infeasible, redistribute the control effect to the remaining free
  % variables, if there are any.
  while any([i_min ; i_max]) & any(i_free) & (iter<imax);
    iter = iter + 1;
    % Now solve for optimal values of the remaining free variables.
    % See 2002-02-07.
    Wu1 = Wu(i_free,:);
    Wu11 = Wu1(:,i_free);
    B1 = B(:,i_free);
    
    A = Wv*B1/Wu11;
    b = Wv*(v-B*u-B1*(Wu11\(Wu1*(ud-u))));
    % Solve for optimal perturbation.
    p1 = Wu11\(pinv_sol(A,b)+Wu1*(ud-u));
    % Update solution
    u(i_free) = u(i_free) + p1;
    
    % Find indeces of _new_ infeasible components.
    i_min = u < umin;
    i_max = u > umax;
    % Set these variables to their limits.
    u(i_min) = umin(i_min);
    u(i_max) = umax(i_max);
    % Remove these from the free variables.
    i_free = i_free & ~(i_min | i_max);
  end
