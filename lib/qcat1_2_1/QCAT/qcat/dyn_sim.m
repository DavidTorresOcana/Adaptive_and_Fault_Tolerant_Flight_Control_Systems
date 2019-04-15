function [u,W,time,iter] = dyn_sim(B,v,plim,varargin)

% DYN_SIM - Dynamic control allocation simulation. 
%
%  [u,W,time,iter] = dyn_sim(B,v,plim,[rlim,T,Wv,W1,W2,S],options)
%
% Performs dynamic control allocation for a sequence of virtual
% control commands stored in v. For each value of v, the control
% signal u is determined by solving
%  
%   min   ||W1(u(t)-Sv(t))||^2 + ||W2(u(t)-u(t-T))||^2
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
% B     control effectiveness matrix (k x m)
% v     commanded virtual control trajectory (k x N)
% plim  position limits [min max] (m x 2)
% rlim  rate limits [min max] (m x 2) ([] --> no rate limits)
% T     sampling time [1]
% Wv    virtual control weighting matrix (k x k) [I]
% W1    control position weighting matrix (m x m) [I]
% W2    control rate weighting matrix (m x m) [0]
% S     steady state control distribution (m x k) [0]
%
%  Options: See QP_SIM
%  --------
%
%  Outputs:
%  -------
% u     optimal control
% W     active constraints in u (+/- : max/min, 1/2 : position/rate)
% time  average computation time per sample
% iter  no. of iterations
%
%  Step response example:
%
%   B = [2 1]; t = 0:.2:10; v = 1*(t>1); plim = [-1 1;-1 1];
%   W1 = eye(2); W2 = diag([5 0]); S = pinv(B);
%   u = dyn_sim(B,v,plim,[],.2,1,W1,W2,S);
%   figure(1),bodemag(dca(B,S,W1,W2,.2))
%   figure(2),stairs(t,[u' v']),legend('u_1','u_2','v=2u_1+u_2')
%
% See also: DCA, QP_SIM, DIR_SIM.
  
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
  W1   = eye(m);
  W2   = zeros(m);
  S    = zeros(m,k);

  % Set values of submitted optional arguments
  for i=1:narg
    switch i
     case 1, rlim = varargin{i};
     case 2, T	  = varargin{i};
     case 3, Wv	  = varargin{i};
     case 4, W1	  = varargin{i};
     case 5, W2	  = varargin{i};
     case 6, S	  = varargin{i};
    end
  end
  
  % Call generic control allocation simulation subroutine.
  [u,W,time,iter] = alloc_sim('dyn',B,v,plim,rlim,T,'Wv',Wv,'W1', ...
			      W1,'W2',W2,'S',S,varargin{iopt:end});
  