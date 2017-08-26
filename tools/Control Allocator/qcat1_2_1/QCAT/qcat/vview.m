function ratio = vview(B,plim,P)
  
% VVIEW - View the attainable virtual control set.
% 
%  1) vview(B,plim)
%
% Shows the attainable virtual control set considering actuator
% position constraints, given by { v : v = B*u, umin < u < umax }.
%
%  2) ratio = vview(B,plim,P)
% 
% Compares the set of feasible virtual control inputs when
%
%  a) the actuator redundancy is fully utilized (as above) [blue]
%  b) a linear allocation control law u = Pv is used (BP = I) [red]
%
% The second set is given by { v : umin < P*v < umax }.
%
%  Inputs:
%  -------
% B      control effectiveness matrix (k x m)
% plim   control position limits [min max] (m x 2)
% P      virtual control law matrix (m x k)
% 
%  Outputs:
%  --------
% ratio  The ratio between the sizes (areas, volumes, ...)
%        of the two sets
%
% The result is only graphically illustrated for k = 1, 2, or 3.
%
% See also: VVIEW_DEMO
  
% Model dimensions
  [k,m] = size(B);
  
  % ------------------------------------------------
  %  a) Find maximum attainable virtual control set 
  %     considering constraints.
  % ------------------------------------------------
  
  % Generate matrix to index corners of feasible control set.
  idx = zeros(2^m,m);
  M = 1:m;
  for i = 1:2^m;
    cbin = dec2bin(i-1,m); % '001'
    c = str2num(cbin')'; % [0 0 1]
    c = c(end:-1:1); % [1 0 0]
    idx(i,:) = 2*M - c;
  end

  % Generate corner points of the feasible control set.
  plimT = plim';
  U = plimT(idx)';

  % Compute the corresponding points in the virtual control space
  V = B*U;

  if nargin > 2
    
    % ---------------------------------------------
    %  b) Find attainable virtual control set when
    %     a linear control law u=Pv is used.
    % ---------------------------------------------
    
    % We want to determine where the k-dim. hyperplane Pv
    % intersects the m-dim. hyperbox of feasible controls.
    % To get the corner points of this set, solve
    % Pv = x where x has k specified entries.
    % 
    % Example: m=3, k=1 -> points will lie on surfaces
    %          m=3, k=2 -> points will lie on edges
    
    % Generate index matrix for all combinations of min and max indeces
    % in k dimensions.
    sub_idx = idx(1:2^k,1:k);
    
    Ulin = [];
    % Loop over all combinations of dimensions
    i_dim = nchoosek(1:m,k);
    for i = 1:size(i_dim,1)
      % For each combination, compute the intersections with all
      % possible min/max combinations.
      
      % k-dimensional min/max combinations
      sub_plimT = plimT(:,i_dim(i,:));
      sub_u_boundary = sub_plimT(sub_idx)';
      
      % Determine which virtual control sub_u_boundary corresponds to
      sub_P = P(i_dim(i,:),:);
      if rank(sub_P) == k % Avoid "parallel" cases
        % Solve sub_u_boundary = sub_P v for v
	v = sub_P\sub_u_boundary;
	% Determine the full countol vector (contains sub_u_boundary)
	u_boundary = P*v;
	
	% Store feasible points
	i_feas = feasible(u_boundary,plim);
	Ulin = [Ulin u_boundary(:,i_feas)];
      end
    end

    % Compute the corresponing points in the virtual control space
    Vlin = B*Ulin;
  
  end
    
  % Compute and visualize the convex hull of the set(s)
  clf
  switch k
   case 1
    K = [min(V) max(V)];
    if nargin > 2
      Klin = [min(Vlin) max(Vlin)];
      ratio = diff(Klin)/diff(K);
      
      % Illustrate
      plot(K,[0 0],'b-o',Klin,-[0 0],'r-o')
    else
      plot(K,[0 0],'b-o')
    end
    xlabel('v')
    
   case 2
    [K,area1] = convhull(V(1,:),V(2,:));
    if nargin > 2
      [Klin,area2] = convhull(Vlin(1,:),Vlin(2,:));
      ratio = area2/area1;
      
      % Illustrate
      fill(V(1,K),V(2,K),[.95 .95 1],...
	   Vlin(1,Klin),Vlin(2,Klin),[1 1 .9])
      hold on;
      plot(Vlin(1,Klin),Vlin(2,Klin),'r',V(1,K),V(2,K),'b')
      hold off;
    else
      fill(V(1,K),V(2,K),[.95 .95 1]);
      hold on;
      plot(V(1,K),V(2,K),'b')
      hold off;
    end
    axis equal;
    xlabel('v_1')
    ylabel('v_2')
    
   otherwise
    [K,vol1]    = convhulln(V');
    if nargin > 2
      [Klin,vol2] = convhulln(Vlin');
      ratio = vol2/vol1;
    end
      
    if k == 3
      % Illustrate
      if nargin > 2
	h = polyplot(Klin,Vlin',1);
	set(h,'EdgeColor','r','FaceColor',[1 1 .9]);
	hold on;
	% Fix: Make V wireframe enclose Vlin
	V0 = repmat(mean(V')',1,size(V,2));
	V = 1.0001*(V-V0)+V0;
	h = polyplot(K,V',1);
	set(h,'EdgeColor',[.6 .6 1],'FaceColor','none');
	hold off
      else
	h = polyplot(K,V',1);
	set(h,'EdgeColor','b','FaceColor',[.95 .95 1]);
      end
      
      xlabel('v_1')
      ylabel('v_2')
      zlabel('v_3')
      view(3);
      axis equal;
      axis vis3d;
      grid on;
    end
  end
  
function f = feasible(x,plim)
% x   m*n
% lb  m
% ub  m
  
  m = size(x,1);
  
  % Mean point
  x0 = mean(plim,2);
  
  % Make the mean point the origin
  x = x - x0*ones(1,size(x,2));
  lb = plim(:,1) - x0; % < 0
  ub = plim(:,2) - x0; % > 0
  
  % Check for feasibility
  tol = 1e-5;
  f = sum((diag(1./ub)*x <= 1+tol) & (diag(1./lb)*x <= 1+tol)) == m;

function h = polyplot(face,vert,merge)
  
  if merge 
    % Merge adjacent, parallel triangles to get fewer lines that
    % are not edges of the polyhedron.
    face4 = [];
    % Loop over all combinations of triangles
    k = 1;
    while k < size(face,1)
      l = k+1;
      while l <= size(face,1)
	iv = intersect(face(k,:),face(l,:)); % Intersecting vertices
	if length(iv) == 2 % Two common vertices
	  % Are the faces parallel?
	  niv = setxor(face(k,:),face(l,:)); % Non-intersecting vertices
	  % Vectors from first common vertex to remaining three vertices
	  A = [vert(iv(2),:)  - vert(iv(1),:);
	       vert(niv(1),:) - vert(iv(1),:);
	       vert(niv(2),:) - vert(iv(1),:)];
	  if abs(det(A))<100*eps
	    % Vectors lie in same plane -> create patch with four vertices
	    face4 = [face4 ; iv(1) niv(1) iv(2) niv(2)];
	    % ... and remove the two triangles
	    face = face([1:k-1 k+1:l-1 l+1:end],:);
	    k = k-1;
	    break
	  end	  
	end
	l = l+1;
      end % inner loop
      k = k+1;
    end % outer loop
    h = [patch('Faces',face,'Vertices',vert)
	 patch('Faces',face4,'Vertices',vert)];
  else
    % Just plot the polyhedron made up by triangles
    h = patch('Faces',face,'Vertices',vert);
  end
  
  