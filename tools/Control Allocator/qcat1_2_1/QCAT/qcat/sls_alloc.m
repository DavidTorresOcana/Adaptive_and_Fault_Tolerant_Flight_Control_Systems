function [u,W,iter] = sls_alloc(B,v,umin,umax,Wv,Wu,ud,u,W,imax)
  
% SLS_ALLOC - Control allocation using sequential least squares.
%
%  [u,W,iter] = sls_alloc(B,v,umin,umax,[Wv,Wu,ud,u0,W0,imax])
%
% Solves the bounded sequential least-squares problem
%
%   min ||Wu(u-ud)||   subj. to   u in M
%
% where M is the set of control signals solving
%
%   min ||Wv(Bu-v)||   subj. to   umin <= u <= umax
%
% using a two stage active set method. To handle the case of coplanar
% controls, pass B+i instead of just B.
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
% See also: MLS_ALLOC, CGI_ALLOC, QP_SIM, ISCOPLANAR.
  
% Numerical tolerance
  tol = 1e-10;
  
  % Check if the coplanar safe version of the algorithm should be
  % used.
  coplanar = ~isreal(B);
  B = real(B);
  
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
  % Compute "A"-matrix and initial residual:
  % ||Wv(B(u0+p)-v)|| = ||Ap-d||
  A = Wv*B;
  d = Wv*(v-B*u);
  
  % Iterate until optimum is found or maximum number of iterations
  % is reached.
  for iter = 1:imax
    % ----------------------------------------
    %  Compute optimal perturbation vector p.
    % ----------------------------------------
    
    % Determine indeces of free variables
    i_free = W==0;
    % and number of free variables.
    m_free = sum(i_free);
    
    if phase == 1
      % ---------
      %  Phase 1 
      % ---------
      
      % Eliminate saturated variables.
      A_free = A(:,i_free);
      
      % Solve the reduced optimization problem for free variables. If
      % A_free p_free = d does not have a unique least-squares
      % solution, pick the minimum length solution.
      if coplanar
	% Actuators are allowed to produce coplanar controls. This means
	% that A_free may not be of full rank even if m_free>k.
	p_free = pinv_sol(A_free,d);
      else
	% A_free will allways be of full (row or column) rank. This
	% leads to two different cases:
	if m_free <= k
	  % A_free p_free = d has a unique least-squares solution.
	  p_free = A_free\d;
	else
	  % Pick minimum norm solution to A_free p_free = d.
	  [Q1,R1] = qr(A_free',0);
	  p_free = Q1*(R1'\d);
	end
      end 
	  
      % Insert perturbations from p_free into free the variables.
      p = zeros(m,1);
      p(i_free) = p_free;
      
    else
      % ---------
      %  Phase 2 
      % ---------

      % Determine indeces of fixed variables,
      i_fixed = find(W);
      % and number of fixed variables.
      m_fixed = length(i_fixed);
      % Construct C0 containing the active inequalities.
      C0 = zeros(m_fixed,m);
      for i = 1:m_fixed
	C0(i,i_fixed(i)) = -W(i_fixed(i));
      end
      % Construct E containing all equality constraints. By
      % construction, the rows of E (= the equality constraints) are
      % linearly independent.
      E = [B ; C0];

      % Number of  constraints.
      k_c = k + m_fixed;
      
      % Compute its QR decomposition E' = QR =(Q1 Q2)(R1)
      %                                              ( 0)
      % Note that the computation of lambda also uses this
      % decomposition.
      [Q,R] = qr(E');
      Q1 = Q(:,1:k_c);
      Q2 = Q(:,k_c+1:end);
      R1 = R(1:k_c,:);
      % Optimal solution:      
      q2 = (A*Q2)\d;
      p = Q2*q2;
      
    end 
    % Optimal perturbation p computed.
    
    % --------------------------------
    %  Is the optimal point feasible?
    % --------------------------------
    
    u_opt = u + p;
    if ~any(u_opt(i_free) < umin(i_free) | u_opt(i_free) > umax(i_free))
      
      % ------
      %  Yes.
      % ------
      
      % Update point.
      u = u_opt;
      % Update residual.
      d = d - A*p;
      
      % --------------------------------------------------------
      %  Is the point the optimal solution to the full problem?
      % --------------------------------------------------------
    
      if (~coplanar) & (phase == 1) & (m_free >= k)
	% No planar controls. If u is feasible, then Bu=v must hold,
        % by construction, i.e., one phase 1 optimum found.
	% Move to phase 2.
	phase = 2;
	% Update "A"-matrix and residual.
	A = Wu;
	d = A*(ud-u);
      else
	% Check for optimality by computing the Lagrange multipliers.
	
	% Compute gradient of optimization criterion.
	g = -A'*d; % d = b-Au --> gradient g = A'(Au-b)
	
	% Compute lambda for all bounds (including inactive ones).
	if phase == 1
	  % Only box constraints.
	  lambda = -W.*g;
	else
	  % Box constraints and equality constraints. Use QR decomp.
	  Lambda = R1\(Q1'*g);
	  % Extract multipliers related to bounds.
	  lambda = zeros(m,1);
	  lambda(i_fixed) = Lambda(k+1:end);
	end
	
	if lambda >= -tol
	  % All relevant Lagrange multipliers are positive.
	  if coplanar & (phase == 1)
	    % Although ||Wv(Bu-v)|| is minimal, ||Wu(u-ud)|| may be
            % reduced further in phase 2.
	    phase = 2;
	    % Update "A"-matrix and residual.
	    % ||Wu(u+p-ud)|| = ||Ap-d||
	    A = Wu;
	    d = Wu*(ud-u);
	    % Caveat: Adding the constraint Bp = 0 in phase 2 may
            % cause linear dependence among the active constraints.
	    % Therefore remove the active box constraints causing
            % such linear dependence.
	    %
	    % Not found a way to do this with quickly with
            % precision --> use "brute force" and empty the working
            % set.
	    W = zeros(m,1);
	  else
	    % / ------------------------ \
	    % | Optimum found, bail out. |
	    % \ ------------------------ /
	    return;
	  end
	else
	  
	  % --------------------------------------------------
	  %  Optimum not found, remove one active constraint.
	  % --------------------------------------------------
	  
	  % Remove constraint with most negative lambda from the
	  % working set.
	  [lambda_neg,i_neg] = min(lambda);
	  W(i_neg) = 0;
	end % lambda >= 0
      end	
    else % feasible?
      
      % ------------------------------------------------------------
      %  No (point not feasible), find primary bounding constraint.
      % ------------------------------------------------------------
      
      % Compute distances to the different boundaries. Since alpha < 1
      % is the maximum step length, initiate with ones.
      dist = ones(m,1);
      i_min = i_free & p<-tol; % rather than <0
      i_max = i_free & p>tol;
      
      dist(i_min) = (umin(i_min) - u(i_min)) ./ p(i_min);
      dist(i_max) = (umax(i_max) - u(i_max)) ./ p(i_max);
      
      % Proportion of p to travel.
      [alpha,i_alpha] = min(dist);
      % Update point and residual.
      u = u + alpha*p;
      d = d - A*alpha*p;
      
      % Add corresponding constraint to working set.
      W(i_alpha) = sign(p(i_alpha));
      
    end
    
  end
