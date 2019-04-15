function [u,W,time,iter] = qp_sim(B,v,plim,varargin)

% QP_SIM - QP control allocation simulation. 
%
%  [u,W,time,iter] = qp_sim(B,v,plim,[rlim,T,Wv,Wu,ud],options)
%
% Performs l2-optimal control allocation for a sequence of virtual
% control commands stored in v. For each value of v, the control
% signal u is determined by solving the quadratic program
%  
%   min ||Wu(u-ud)||   subj. to   u in M
%
% where M is the set of control signals solving
%
%   min ||Wv(Bu-v)||   subj. to   u in U
%
% where U is the set of feasible control signals with respect to
% position limits and, optionally, rate limits.
%
%  Inputs:
%  -------
% B     control effectiveness matrix (k x m)
% v     commanded virtual control trajectory (k x N)
% plim  position limits [min max] (m x 2)
% rlim  rate limits [min max] (m x 2) ([] --> no rate limits)
% T     sampling time [1]
% Wv    virtual control weighting matrix (k x k) [I]
% Wu    control weighting matrix (m x m) [I]
% ud    desired control (m x 1) [0]
%
%  Options: options = option_1,value_1,option_2,value_2,...
%  --------
% 'alg'    numerical algorithm: 'sls'    SLS_ALLOC
%                               'wls'    WLS_ALLOC (default)
%                               'wlsc'   WLSC_ALLOC
%                               'mls'    MLS_ALLOC
%                               'ip'     IP_ALLOC
%                               'cgi'    CGI_ALLOC
%                               'fxp'    FXP_ALLOC
% 'imax'   max no. of iterations [100]
% 'gam'    weight used in algorithms based on weighted LS [1e6]
% 'tol'    tolerance used in IP_ALLOC stopping criterion [1e-6]
% 'hot'    hotstart solver (not ip/cgi) with previous solution (0/[1])
% 'ui'     initial control signal
% 'Wi'     initial active constraints
% 'rep'    no. of repetitions [1]
%
%  Outputs:
%  -------
% u     optimal control
% W     active constraints in u (+/- : max/min, 1/2 : position/rate)
% time  average computation time per sample
% iter  no. of iterations
%
% Flight example:
%   load admiredata
%   u=qp_sim(B,v,plim,rlim,T);
%   figure(1),plot(t,u*180/pi),legend('can','elev (r)','elev (l)','rud')
%   figure(2),plot(t,B*u,t,v,'k--'),legend('roll','pitch','yaw')
%
% See also: DIR_SIM, DYN_SIM, and the allocation algorithms above.

% F18 example
%   load f18data
%   u=qp_sim(B,v,plim,rlim,T1);
%   plot(tn,u)
%   plot(tn,v,'k',tn,B*u)
  
% Number of variables
  [k,m] = size(B);

  % Find no. of optional arguments (excluding options)
  iopt = length(varargin)+1;
  for i = 1:length(varargin)
    if ischar(varargin{i})
      iopt = i; % index of first option string
      break;
    end
  end
  narg = iopt-1;
  
  % Set default values of optional arguments
  rlim = [];
  T    = 1;
  Wv   = eye(k);
  Wu   = eye(m);
  ud   = zeros(m,1);

  % Set values of submitted optional arguments
  for i=1:narg
    switch i
     case 1, rlim = varargin{i};
     case 2, T	  = varargin{i};
     case 3, Wv	  = varargin{i};
     case 4, Wu	  = varargin{i};
     case 5, ud	  = varargin{i};
    end
  end

  % Call generic control allocation simulation subroutine.
  [u,W,time,iter] = alloc_sim('qp',B,v,plim,rlim,T,'Wv',Wv,...
			      'Wu',Wu,'ud',ud,varargin{iopt:end});
  