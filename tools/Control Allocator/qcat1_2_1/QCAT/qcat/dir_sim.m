function [u,W,time] = dir_sim(B,v,plim,varargin)

% DIR_SIM - Direct control allocation simulation. 
%
%  [u,W,time] = dir_sim(B,v,plim,[rlim,T],options)
%
% Performs direct control allocation for a sequence of virtual
% control commands stored in v. For each value of v, the control
% signal u is determined by solving
%
%   max a   subj. to  Bu = av
%   a,u               u in U
%
% If a > 1,  u = u/a.
%
% Here, U is the set of feasible controls only with respect to
% the position limits. Rate limiting is (optionally) performed
% on the allocated control vector.
%
%  Inputs:
%  -------
% B     control effectiveness matrix (k x m)
% v     commanded virtual control trajectory (k x N)
% plim  position limits [min max] (m x 2)
% rlim  rate limits [min max] (m x 2) ([] --> no rate limits)
% T     sampling time [1]
%
%  Options: options = option_1,value_1,option_2,value_2,...
%  --------
% 'ui'  initial control signal
% 'rep' no. of repetitions [1]
%
%  Outputs:
%  -------
% u     optimal control
% W     active constraints in u (+/- : max/min, 1/2 : position/rate)
% time  average computation time per sample
% 
%  Example:
%
%   load f18data
%   u=dir_sim(B,v,plim,rlim,T2);
%   figure(1),plot(tn,u*180/pi),ylabel('Controls (deg)')
%   figure(2),plot(tn,B*u,tn,v,'k--'),legend('roll','pitch','yaw')
%
% See also: DIR_ALLOC, QP_SIM, DYN_SIM.
  
% Different example:
%   B = [2 1]; t = 0:.2:10; v = 4*sin(t); plim = [-1 1;-2 1];
%   u = dir_sim(B,v,plim);
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

  % Set values of submitted optional arguments
  for i=1:narg
    switch i
     case 1, rlim = varargin{i};
     case 2, T	  = varargin{i};
    end
  end
      
  % Call generic control allocation simulation subroutine.
  [u,W,time] = alloc_sim('dir',B,v,plim,rlim,T,'hot',0, ...
			 varargin{iopt:end});
  