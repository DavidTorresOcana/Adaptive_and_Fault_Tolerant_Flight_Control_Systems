function [c,cols] = iscoplanar(B)
  
% ISCOPLANAR - Test for coplanar controls.
%
%  [c,cols] = iscoplanar(B)
% 
% Given a control effectiveness matrix B (k x m), the function returns
% c = 1 if there are k columns in B that are linearly dependent.
% The columns indexes of all such combinations are gathered in cols.
% 
% Example:
%
%   B = [1 1 0;0 0 1]
%   [c,cols] = iscoplanar(B)
  
% Dimensions
  [k,m] = size(B);
  
  % Construct all n combinations of k columns.
  C = nchoosek(1:m,k);
  n = nchoosek(m,k);
  
  c = 0;
  cols = [];
  % Check all k x k submatrices of B.
  for i = 1:n
    Bsub = B(:,C(i,:));
    if qrank(Bsub) < k
      % Coplanar controls detected
      c = 1;
      cols = [cols ; C(i,:)];
    end
  end