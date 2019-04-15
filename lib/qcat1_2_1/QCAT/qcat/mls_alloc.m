function [u,W,iter] = mls_alloc(B,v,umin,umax,Wv,Wu,ud,u,W,imax)
  
% MLS_ALLOC - Control allocation using minimal least squares.
%
%  [u,W,iter] = mls_alloc(B,v,umin,umax,[Wv,Wu,ud,u0,W0,imax])
%
% Solves the bounded sequential least-squares problem
%
%   min ||Wu(u-ud)||   subj. to   u in M
%
% where M is the set of control signals solving
%
%   min ||Wv(Bu-v)||   subj. to   umin <= u <= umax
%
% using a two stage active set method. Wu must be diagonal since the
% problem is reformulated as a minimal least squares problem. The
% implementation does not handle the case of coplanar controls.
%
%  Inputs:
%  -------
% B     control effectiveness matrix (k x m)
% v     commanded virtual control (k x 1)
% umin  lower position limits (m x 1)
% umax  upper position limits (m x 1)
% Wv    virtual control weighting matrix (k x k) [I]
% Wu    control weighting matrix (m x m), diagonal [I]
% ud    desired control (m x 1) [0]
% u0    initial point (m x 1)
% W0    initial working set (m x 1) [empty]
% imax  max no. of iterations [100]
% 
%  Outputs:
%  -------
% u     optimal control
% W     optimal active set
% iter  no. of iterations (= no. of changes in the working set + 1)
%
%                           0 if u_i not saturated
% Active set syntax: W_i = -1 if u_i = umin_i
%                          +1 if u_i = umax_i
%
% See also: SLS_ALLOC, CGI_ALLOC, QP_SIM.
  
% k = number of virtual controls
% m = number of variables (actuators)
  [k,m] = size(B);
  
  % Set default values of optional arguments
  if nargin < 10
    imax = 100; % Heuristic value
    if nargin < 9, u = (umin+umax)/2; W = zeros(m,1); end
    if nargin < 7,  ud = zeros(m,1); end
    if nargin < 6,  Wu = eye(m);     end
    if nargin < 5,  Wv = eye(k);     end
  end

  % Start with phase one.
  phase = 1;
  
  % Reformulate as a minimal least squares problem. See 2002-03-08 (1).
  A = Wv*B/Wu;
  b = Wv*(v-B*ud); % Note sign convention: min ||Ax-b||
  xmin = Wu*(umin-ud);
  xmax = Wu*(umax-ud);
  % Compute initial point and residual.
  x = Wu*(u-ud);
  r = A*x-b;
  % Determine indeces of free variables
  i_free = W==0;
  % and number of free variables.
  m_free = sum(i_free);
  
  % Iterate until optimum is found or maximum number of iterations
  % is reached.
  for iter = 1:imax
    % ----------------------------------------
    %  Compute optimal perturbation vector p.
    % ----------------------------------------
    
    if phase == 1
      % ---------
      %  Phase 1 
      % ---------

      % Eliminate saturated variables.
      A_free = A(:,i_free);

      % Solve the reduced optimization problem for free variables.
      % Under the assumption that no n actuators produce coplanar
      % control efforts, A_free will allways be of full rank. This
      % leads to two different cases:
      if m_free <= k
	% --------------------------------------------------
	%  min ||A_free p_free + r|| has a unique solution.
	% --------------------------------------------------
	p_free = -A_free\r;
      else
	% ---------------------------------------------------
	%  min ||A_free p_free + r|| has no unique solution.
	%
	%  Pick the minimal solution.
	% ---------------------------------------------------

	% QR decompose: A_free' = QR = (Q1 Q2)(R1)
	%                                     ( 0)
	[Q1,R1] = qr(A_free',0);
	p_free = -Q1*(R1'\r);
      end
      
      % Insert perturbations from p_free into free the variables.
      p = zeros(m,1);
      p(i_free) = p_free;

    else
      % ---------
      %  Phase 2 
      % ---------
      
      % Determine indeces of fixed variables,
      i_fixed = ~i_free;
      % and number of fixed variables.
      m_fixed = m - m_free;

      % HT' = rows of U in working set.
      % HT = (m - n) x i_fixed
      HT = U(i_fixed,:)';
      % QR decompose: HT = VRtot = (V1 V2)(R)
      %				          (0)
      % Note that the computation of lambda also uses this
      % decomposition.
      [V,Rtot] = qr(HT);
      V1 = V(:,1:m_fixed);
      V2 = V(:,m_fixed+1:end);
      R = Rtot(1:m_fixed,:);
      % Optimal solution:
      s = -V2'*z;
      pz = V2*s;
      p = U*pz;
      
    end % Optimal perturbation p computed.
    
    % --------------------------------
    %  Is the optimal point feasible?
    % --------------------------------
    
    x_opt = x + p;
    infeasible = (x_opt < xmin) | (x_opt > xmax);

    if ~any(infeasible(i_free))
    
      % ------
      %  Yes.
      % ------
      
      % Update point and residual
      x = x_opt;
      if phase == 1
	r = r + A*p;
      else
	z = z + pz;
      end
      
      if (phase == 1) & (m_free >= k)
	% If u is feasible, then Bu=v must hold, by construction.
	% Move to phase 2.
	phase = 2;

	% QR decompose A'.
	[Utot,Stot] = qr(A');
	U = Utot(:,k+1:end);
	% Initial point.
	z = U'*x;
	
      else
	
	% Check for optimality by computing the Lagrange multipliers.
	% Compute lambda for all bounds (including inactive ones).
	lambda = zeros(m,1);

	if m_free < m
	  if phase == 1
	    % Compute gradient of optimization criterion.
	    g = A'*r;
	    % Compute Lagrange variables.
	    lambda = -W.*g;
	    
	  else
	    % Insert multipliers related to active bounds.
	    lambda(i_fixed) = -W(i_fixed).*(R\(V1'*z));
	  end
	end
	  
	if lambda >= -eps
	  % / ------------------------ \
	  % | Optimum found, bail out. |
	  % \ ------------------------ /

	  % Compute solution in original coordinates.
	  u = Wu\x + ud;
	  return;
	end
	
	% --------------------------------------------------
	%  Optimum not found, remove one active constraint.
	% --------------------------------------------------
	
	% Remove constraint with most negative lambda from the
        % working set.
	[lambda_neg,i_neg] = min(lambda);
	W(i_neg) = 0;
	i_free(i_neg) = 1;
	m_free = m_free + 1;
      end
    else
      
      % ------------------------------------------------------------
      %  No (point not feasible), find primary bounding constraint.
      % ------------------------------------------------------------
      
      % Compute distances to the different boundaries. Since alpha < 1
      % is the maximum step length, initiate with ones.
      dist = ones(m,1);
      i_min = i_free & p<0;
      i_max = i_free & p>0;
      
      dist(i_min) = (xmin(i_min) - x(i_min)) ./ p(i_min);
      dist(i_max) = (xmax(i_max) - x(i_max)) ./ p(i_max);
      
      % Proportion of p to travel.
      [alpha,i_alpha] = min(dist);

      % Update point (and residual).
      x = x + alpha*p;
      if phase == 1
	r = r + A*alpha*p;
      else
	z = z + alpha*pz;
      end
	
      % Add corresponding constraint to working set.
      W(i_alpha) = sign(p(i_alpha));
      i_free(i_alpha) = 0;
      m_free = m_free - 1;
      
    end
    
  end
  
  % Compute solution in original coordinates.
  u = Wu\x + ud;