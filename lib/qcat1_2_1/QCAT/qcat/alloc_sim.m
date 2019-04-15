function [u,A,time,iter] = alloc_sim(method,B,v,plim,rlim,T,varargin)

% ALLOC_SIM - Control allocation simulation.
%
%  [u,A,iter,time] = alloc_sim(method,B,v,plim,rlim,T,varargin)
%
% Subroutine not intended to be called directly by user.
  
% Performs control allocation for a sequence of virtual control
% commands stored in v.
%
%  Inputs:
%  -------
% method   control allocation method: 'qp'   l2-optimal allocation
%				      'dyn'  dynamic allocation
%				      'dir'  direct allocation
% B	   control effectiveness matrix (k x m)
% v        commanded virtual control trajectory (k x N)
% plim     position limits [min max] (m x 2)
% rlim	   rate limit = [min max] (m x 2) 
% T	   sampling time
% options  options, see code
%
%  Outputs:
%  -------
% u     optimal control
% A     active constraints in u
% iter  no. of iterations
% time  average computation time per sample
  
% Number of variables
  [k,m] = size(B);
  % Coplanar controls? Only in SLS_ALLOC
  copl = iscoplanar(B);
  % Number of samples
  N = size(v,2);
  % Allocate output variables
  u = zeros(m,N);
  A = zeros(m,N);
  iter = zeros(1,N);
  time = zeros(1,N);
  
  % Set default values of optional arguments
  switch method
   case {'qp','dyn'}
    alg  = 'wls';
   case 'dir'
    alg = 'dir';
  end
  M    = 1;	     % no of repetitions
  imax = 100;	     % no of iterations
  gam  = 1e6;	     % weight
  tol  = 1e-6;       % tolerance in IP solver
  hs   = 1;	     % hotstart
  ui   = [];	     % initial control
  Wi   = zeros(m,1); % initial working set
  Wv   = eye(k);     % QP allocation
  Wu   = eye(m);
  ud   = zeros(m);
  W1   = eye(m);     % Dynamic allocation
  W2   = zeros(m);
  S    = pinv(B);
  
  for i=1:2:length(varargin)
    switch(varargin{i})
      % --- General options ---
     case 'alg'  % solver, algorithm
      alg = varargin{i+1};
     case 'rep'  % no of repetitions
      M = varargin{i+1};

      % --- Method specific options ---
     case 'imax' % number of iterations
      imax = varargin{i+1};
     case 'gam'  % LS weight
      gam = varargin{i+1};
     case 'tol'  % IP tolerance
      tol = varargin{i+1};
     case 'hot'  % hotstart
      hs = varargin{i+1};
     case 'ui'   % active set (also fix-point)
      ui = varargin{i+1};
     case 'Wi'
      Wi = varargin{i+1};
     case 'Wv'   % QP allocation
      Wv = varargin{i+1};
     case 'Wu'
      Wu = varargin{i+1};
     case 'ud'
      ud = varargin{i+1};
     case 'W1'   % dynamic allocation
      W1 = varargin{i+1};
     case 'W2'
      W2 = varargin{i+1};
     case 'S'
      S = varargin{i+1};
     
     otherwise
      error(sprintf('Unknown option: %s',varargin{i}));
    end
  end

  if strcmp(alg,'ip')
    % Weighting matrices must be unit matrices.
    if any(any(Wv ~= eye(k))) | any(any(Wu ~= eye(m)))
      disp(' ')
      disp(['** Warning: Non-unit matrices Wv and Wu not handled by' ...
	    ' IP_ALLOC **']) 
      disp(' ')
    end
    % Dynamic allocation almost always requires non-unit matrix Wu.
    if strcmp(method,'dyn')
      disp(' ')
      disp(['** Warning: Dynamic allocation not handled by IP_ALLOC,' ...
	    ' resorting to WLS_ALLOC **'])
      disp(' ');
      alg = 'wls';
    end
  end
  
  if isempty(ui)
    % Set default initial conditions
    switch method
     case 'qp'
      [ui,Wi] = wls_alloc(B,v(:,1),plim(:,1),plim(:,2),Wv,Wu,ud);
     case 'dyn'
      [ui,Wi] = wls_alloc(B,v(:,1),plim(:,1),plim(:,2),Wv,W1,S*v(:,1));
     case 'dir'
      ui = dir_alloc(B,v(:,1),plim(:,1),plim(:,2));
    end
  end
  
  % Hotstart?
  if hs
    u0 = ui;
    W0 = Wi;
  end
  
  % Precompute matrices used in dynamic allocation
  if strcmp(method,'dyn')
    W1sq = W1^2;
    W2sq = W2^2;
    Wu = sqrtm(W1sq+W2sq);
    invWusq = inv(W1sq+W2sq);
  end
  
  % Allocate control signals to produce the demanded virtual
  % control trajectory.
  uprev = ui;
  for i=1:N
    % Compute feasible upper and lower bounds.
    if isempty(rlim)
      umin = plim(:,1);
      umax = plim(:,2);
    else
      umin = max(plim(:,1),uprev+rlim(:,1)*T);
      umax = min(plim(:,2),uprev+rlim(:,2)*T);
    end
    if ~hs
      u0 = (umin+umax)/2;
      W0 = zeros(m,1);
    end
    % Update u0 to reflect working set. Crucial when rate limits
    % are included and hotstart is used.
    i_min = W0 == -1;
    i_max = W0 == +1;
    u0(i_min) = umin(i_min);
    u0(i_max) = umax(i_max);
    % For dynamic allocation, merge the position and rate terms
    % into one term ||Wu(u-ud)||.
    if strcmp(method,'dyn')
      us = S*v(:,i);
      ud = invWusq*(W1sq*us+W2sq*uprev);
    end
    % Solve control allocation problem M times
    tic;
    for l = 1:M
      switch alg
       case 'sls'
	[u(:,i),W,iter(i)] = ...
	    sls_alloc(B+copl*j,v(:,i),umin,umax,Wv,Wu,ud,u0,W0,imax);
       case 'mls'
	[u(:,i),W,iter(i)] = ...
	    mls_alloc(B,v(:,i),umin,umax,Wv,Wu,ud,u0,W0,imax);
       case 'wls'
	[u(:,i),W,iter(i)] = ...
	    wls_alloc(B,v(:,i),umin,umax,Wv,Wu,ud,gam,u0,W0,imax);
	% --- FoT25 allocation routine, not in standard QCAT ---
       case 'wls_c'
	[u(:,i),W,iter(i)] = ...
	    wls_alloc_c(B,v(:,i),umin,umax,Wv,Wu,ud,gam,u0,W0,imax);
	% ------------------------------------------------------
       case 'wlsc'
	[u(:,i),W,iter(i)] = ...
	    wlsc_alloc(B,v(:,i),umin,umax,Wv,Wu,ud,gam,u0,W0,imax);
       case 'ip'
	[u(:,i),iter(i)] = ...
	    ip_alloc(B,v(:,i),umin,umax,ud,gam,tol,imax);
       case 'fxp'
	u(:,i) = fxp_alloc(B,v(:,i),umin,umax,Wv,Wu,ud,gam,u0,imax);
       case 'cgi'
	u(:,i) = cgi_alloc(B,v(:,i),umin,umax,Wv,Wu,ud,imax);
       case 'dir'
	u(:,i) = dir_alloc(B,v(:,i),plim(:,1),plim(:,2));
       otherwise
	error(sprintf('Unknown allocation algorithm: %s',alg));
      end
    end
    % Register elapsed time
    time(i) = toc/M;
    % Limit control (only necessary with rate limits and direct alloc)
    u(:,i) = max(umin,min(umax,u(:,i)));
    % Determine active constraints in the final point.
    % +/- : max or min
    % 1/2 : position or rate limit
    A(:,i) = - 2*(u(:,i)==umin) + (u(:,i)==plim(:,1)) ...
	     + 2*(u(:,i)==umax) - (u(:,i)==plim(:,2));
    % Update uprev
    uprev = u(:,i);
    % Hotstart?
    if hs
      u0 = u(:,i);
      if any(strcmp(alg,{'sls','wls','wls_c','wlsc','mls'}))
	W0 = W;
      end
    end
  end