function [u,iter] = ip_alloc(B,v,umin,umax,ud,gam,tol,imax)
  
% IP_ALLOC - Control allocation using interior point method.
%
%  [u,iter] = ip_alloc(B,v,umin,umax,[ud,gamma,tol,imax])
%
% Solves the weighted, bounded least-squares problem
%
%   min ||u-ud||^2 + gamma ||Bu-v||^2   (unit weighting matrices)
%
%   subj. to  umin <= u <= umax
%
% using a primal dual interior point method.
%
%  Inputs:
%  -------
% B     control effectiveness matrix (k x m)
% v     commanded virtual control (k x 1)
% umin  lower position limits (m x 1)
% umax  upper position limits (m x 1)
% ud    desired control (m x 1) [0]
% gamma weight (scalar) [1e4]
% tol   tolerance used in stopping criterion [1e-4]
% imax  max no. of iterations [100]
% 
%  Outputs:
%  -------
% u     optimal control
% iter  no. of iterations
%
% See also: WLS_ALLOC, WLSC_ALLOC, FXP_ALLOC, QP_SIM.
%  
% Contributed by John Petersen.

% Set default values of optional arguments
  if nargin < 8
    imax = 100; % Heuristic value
    [k,m] = size(B);
    if nargin < 7, tol = 1e-4; end
    if nargin < 6, gam = 1e4; end
    if nargin < 5, ud = zeros(m,1); end
  end

  % Reformulate   min  ||u-ud||^2 + gamma ||Bu-v||^2  
  %               s.t. umin <= u <= umax
  %
  % as            min  ||Ax-b||^2 + h ||x-xd||^2
  %               s.t. 0 <= x <= xmax
  %
  % where x=u-umin, h=1/gamma, A=B, b=v-B*umin, xd=ud-umin, xmax=umax-umin 
  h = 1/gam; A = B; b = v - B*umin; xd = ud - umin; xmax = umax - umin;
  
  % ||Ax-b||^2 + h ||x-xd||^2 = 1/2x'Hx + c'x + f(xd)
  c = -2*(b'*A + h*xd');
  
  % Solve QP problem.
  [x,iter] = pdq(A,b,c',xmax,h,tol,imax);
  
  % Optimal control.
  u = x + umin;
  
function [x,iter] = pdq(A,b,c,u,wc,tol,imax);
  % Primal dual IP solver
  
  [k,m] = size(A);	% k = #constraints , m = #variables
  As2=A*sqrt(2);
  [s,w,x,z] = startpt(A,b,c,u,wc);
  rho = .9995; sig = 0.1; m2 = 2*m;
  xs = x.*s; wz = w.*z;
  mu = sig*(sum(xs + wz))/m2; % eq. (7)
  nxl=norm(x,1); iHd = 0; 
  ru = 0; rc = 0; rb = 0;
  
  for iter = 1:imax+1
    if (mu < tol) 
      % Close enough to optimum, bail out.
      break;
    end;   
    rxs = (xs - mu);
    iw = 1./w; rwz = (wz - mu);
    ix = 1./x; ixs = ix.*s;
    iwz = iw.*z;
    d = 2*wc + ixs + iwz; 						
    [ds,dw,dx,dz] = direct(d,As2,rb,rc,ru,rxs,rwz,ix,iw,ixs,iwz,z);
    alpha = stepsize(dx,ds,dw,dz,x,s,w,z);
    ralpha = rho*alpha;
    s = s + ralpha*ds;
    w = w + ralpha*dw;
    x = x + ralpha*dx;
    z = z + ralpha*dz;

    xs = x.*s; wz = w.*z;
    gap = sum(xs + wz)/m2;
    mu = min(.1,100*gap)*gap;
  
  end
  
  iter=iter-1;  %True number of iterations is one less

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initial starting point  
function [s,w,x,z] = startpt(A,b,c,u,wc);
  
  m = size(A,2);
  en = ones(m,1);
  % Start at the center of the constraint box
  x = u/2; w = u - x; z = en; s = en; 
  if 0 
    AA = csm(A');  % csm efficiently computes cb'*cb
    for i=1:m
      AA(i,i) = AA(i,i) + wc;
    end % H = 2(A'A+wcI);
  else
    AA = A'*A + wc*eye(m);
  end
  ec = c + 2*AA*x;  %initial residual error used to initialize z,s;
  g = .1;  sg = 1+g;
  i = find(ec>0); s(i) =  sg*ec(i); z(i) =  g*ec(i); % Hx + c + z - s = 0;
  i = find(ec<0); z(i) = -sg*ec(i); s(i) = -g*ec(i);
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute step direction
function [ds,dw,dx,dz] = direct(d,As2,rb,rc,ru,rxs,rwz,ix,iw,ixs,iwz,z);
  
  ixrxs = ix.*rxs; iwrwz = iw.*rwz;
  rr = iwrwz - ixrxs;

  if 0
    iHd = smwf(d,As2);  % smwf is a fast smw for the specific form used here
    dx = iHd*rr;
  else
    % iHd = inv(diag(d)+As2'*As2);
    dx = (diag(d)+As2'*As2)\rr;
  end
  ds = -ixs.*dx - ixrxs;	dw = -dx;
  dz = -iwz.*dw - iwrwz;
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute stepsize
function alpha = stepsize(dx,ds,dw,dz,x,s,w,z);
  i = find(ds<0); as = min(-s(i)./ds(i));	
  i = find(dw<0); aw = min(-w(i)./dw(i));	
  i = find(dx<0); ax = min(-x(i)./dx(i));	
  i = find(dz<0); az = min(-z(i)./dz(i));	
  alpha = min([aw ax as az 1]);	

  
  
  
  
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  csm.m  compute symmetric matrix from X = B*B'
function X = csm(B);
  
  [m,n] = size(B);
  for i=1:m-1
    for j=i+1:m
      X(i,j) = B(i,:)*B(j,:)';
      X(j,i) = X(i,j);
    end;	
  end;
  for i = 1:m;	X(i,i) = B(i,:)*B(i,:)'; end;
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% smwf.m  computes the inverse of (A + B'B) using Sherman-Morrison-Woodbury formula
%    Fast version for specific matrices.
%   		where A is a diagonal matrix, but designated only by a vector
%													i.e. the input is the vector of the diagonal
%				B is non-square matrices 
%   H = inv(A) - inv(A)B'*inv(B*inv(A)B' + I)*B*inv(A);
%		
function		 H = smwf(A,B);
  
  [k,m] = size(B);
  iA = 1./A; iAB = zeros(m,k); BiA = zeros(k,m);  BAB = zeros(k);
  for i = 1:k    iAB(:,i) = iA.*B(i,:)'; end;
  for i = 1:k		BiA(i,:) = B(i,:).*iA';	end;
  for i=1:k-1
    for j=i+1:k
      BAB(i,j) = B(i,:)*iAB(:,j);
      BAB(j,i) = BAB(i,j);
    end;
  end;
  for i = 1:k;	BAB(i,i) = B(i,:)*iAB(:,i); end;
  Q = BAB;
  for i = 1:k;	Q(i,i) = Q(i,i) + 1; end
  QBA = Q\iAB';
  for i=1:m-1
    for j=i+1:m
      ABQBA(i,j) = iAB(i,:)*QBA(:,j);
      ABQBA(j,i) = ABQBA(i,j);
    end;
  end;
  for i = 1:m;  ABQBA(i,i) = iAB(i,:)*QBA(:,i); end
  H = -ABQBA;
  for i = 1:m
    H(i,i) = iA(i) + H(i,i);
  end
  