function x = pinv_sol(A,b,method)
  
% PINV_SOL - Compute pseudoinverse solution.
%
%  x = pinv_sol(A,b,method)
%
% Returns the pseudoinvers solution (= the minimum norm solution) to
% the linear equation Ax = b.

  if nargin < 3
    method = 1;
  end
  
  switch(method)
   case 1
    % -----------------------
    %  Use QR decomposition.
    % -----------------------
    if isempty(A)
      x = [];
      return
    end
  
    [Q1,R1,E] = qr(A',0);
    % Account for rank deficiency. The tolerance is 10 x tolerance
    % used in \. See qr documentation.
    tol = max(size(A))*10*eps*abs(R1(1,1));
    % Compute numerical rank.
    r = sum(abs(diag(R1))>tol);
    % Compute pseudoinverse solution.
    x = Q1(:,1:r)*(R1(1:r,:)'\b(E));
   case 2
    % -------------------------------------
    %  Use MATLAB's pseudoinverse command.
    % -------------------------------------
    x = pinv(A) * b;
   case 3
    % ------------------------
    %  Use SVD decomposition.
    % ------------------------
    % Adopted from pinv.m.
    [n,m] = size(A);
    if n > m
      % Ax = b is overdetermined.
      [U,S,V] = svd(A,0);
    else
      % Ax = b is underdetermined.
      [V,S,U] = svd(A',0);
    end
    s = diag(S);
    r = sum(s > .0000001);
    Sr_pinv = diag(1./s(1:r));
    x = V(:,1:r)*Sr_pinv*U(:,1:r)'*b;
  end