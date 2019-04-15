function r = qrank(A)
  
% QRANK - Quick matrix rank estimation.
%
%  r = qrank(A)
%
% Provides an estimate of rank(A) based on the QR decomposition of
% A. Faster than Matlab's RANK which relies on the SVD.
  
% Compute QR decomposition.
  [Q1,R1,E] = qr(A',0);
  % The tolerance is 10 x tolerance used in \. See qr documentation.
  tol = max(size(A))*10*eps*abs(R1(1,1));
  % Compute numerical rank.
  r = sum(abs(diag(R1))>tol);
