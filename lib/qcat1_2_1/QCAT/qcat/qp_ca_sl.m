function u = qp_ca_sl(arg,B,plim,rlim,T,Wv,Wu,ud,alg,imax,gam,tol)
  
% Wrapper used in the QP control allocation Simulink block.
  
% Dimensions
  [k,m] = size(B);
  
  % Extract nonconstant input arguments
  v = arg(1:k);
  uprev = arg(k+1:end);
  
  % Overall position limits
  if isempty(rlim)
    umin = plim(:,1);
    umax = plim(:,2);
  else
    umin = max(plim(:,1),uprev+rlim(:,1)*T);
    umax = min(plim(:,2),uprev+rlim(:,2)*T);
  end
  u0 = mean([umin umax]')';
  W0 = zeros(m,1);
  
  switch lower(alg)
   case 'sls'
    u = sls_alloc(B+j,v,umin,umax,Wv,Wu,ud,u0,W0,imax);
   case 'mls'
    u =	mls_alloc(B,v,umin,umax,Wv,Wu,ud,u0,W0,imax);
   case 'wls'
    u = wls_alloc(B,v,umin,umax,Wv,Wu,ud,gam,u0,W0,imax);
   case 'wlsc'
    u = wlsc_alloc(B,v,umin,umax,Wv,Wu,ud,gam,u0,W0,imax);
   case 'ip'
    u = ip_alloc(B,v,umin,umax,ud,gam,tol,imax);
   case 'fxp'
    u = fxp_alloc(B,v,umin,umax,Wv,Wu,ud,gam,u0,imax);
   case 'cgi'
    u = cgi_alloc(B,v,umin,umax,Wv,Wu,ud);
  end
