function u = dyn_ca_sl(arg,B,plim,rlim,T,Wv,W1,W2,S)
  
% Wrapper used in the dynamic control allocation Simulink block.
  
% Dimensions
  [k,m] = size(B);
  
  % Extract nonconstant input arguments
  v = arg(1:k);
  uprev = arg(k+1:end);
  
  % Formulate as a regular QP problem: min ||Wu(u-ud)||
  W1sq = W1^2;
  W2sq = W2^2;
  Wu = sqrtm(W1sq+W2sq);
  invWusq = inv(W1sq+W2sq);

  us = S*v;
  ud = invWusq*(W1sq*us+W2sq*uprev);
  
  % Overall position limits
  if isempty(rlim)
    umin = plim(:,1);
    umax = plim(:,2);
  else
    umin = max(plim(:,1),uprev+rlim(:,1)*T);
    umax = min(plim(:,2),uprev+rlim(:,2)*T);
  end
  
  % Use WLS -- fast and robust.
  u = wls_alloc(B,v,umin,umax,Wv,Wu,ud);