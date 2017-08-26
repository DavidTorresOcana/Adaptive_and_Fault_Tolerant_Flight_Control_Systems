function vview_demo
  
% VVIEW_DEMO - Demo of the VVIEW command.
%
% Choose between a number of examples illustrating VVIEW. Issues
% the command VVIEW(B,plim,pinv(B)) for various choices of B.
  
  while 1
  disp(' ')
  disp(' 1. Random control effectiveness matrix B (1x5)')
  disp(' 2.                                       (2x5)')
  disp(' 3.                                       (3x5)')
  disp(' 4. Admire flight data')
  disp(' 5. HARV flight data')
  disp(' ')
  ch = input('Select example (leave blank to exit): ','s');
  disp(' ')
  
  switch ch
    
   case {'1','2','3'} % Randomize
    k = str2num(ch);
    m = 5; % No. of actuators
    B=randn(k,m);
    plim=[-ones(m,1),ones(m,1)];
    ratio = vview(B,plim,pinv(B));
    title('Random B')
    
   case '4' % Admire data - Härkegård PhD thesis, Sec 10.4
    B = [0    -4.24  4.24  1.48;
	 1.65 -1.27 -1.27  0;
	 0    -0.28  0.28 -0.88];
    
    plim = [-55 -30 -30 -30;
	    25  30  30  30]'*pi/180;
    
    ratio = vview(B,plim,pinv(B));
    xlabel('Cl')
    ylabel('Cm')
    zlabel('Cn')
    title('Admire (3000 m, Mach 0.22)')
    
   case '5' % HARV data - Bordignon PhD thesis, Example 5-2
    B = [-4.382e-2 -5.330e-1 -1.100e-2;
	 4.382e-2 -5.330e-1  1.100e-2;
	 -5.841e-2 -6.486e-2  3.911e-3;
	 5.841e-2 -6.486e-2 -3.911e-3;
	 1.674e-2     0     -7.428e-2;
	 -6.280e-2  6.234e-2     0    ;
	 6.280e-2  6.234e-2     0    ;
	 2.920e-2  1.000e-5  3.000e-4;
	 1.000e-5  3.553e-1  1.000e-5;
	 1.000e-2  1.000e-5  1.485e-1]';
    
    plim = [-4.189 1.833;
	    -4.189 1.833;
	    -5.236 5.236;
	    -5.236 5.236;
	    -5.236 5.236;
	    -1.396 7.854;
	    -1.396 7.854;
	    -5.236 5.236;
	    -5.236 5.236;
	    -5.236 5.236]*1e-1;
    ratio = vview(B,plim,pinv(B)); % Ok, same result (13.7%)
    xlabel('Cl')
    ylabel('Cm')
    zlabel('Cn')
    title('HARV (10000 ft, Mach 0.3, \alpha=12.5^o)')
  
    otherwise break
  end  
  
  disp('---------------------------------------------------------------')
  disp(' ');
  disp(' Control effectiveness matrix: B ='),disp(' '),disp(B)
  disp(' Position limits: [umin umax]'' ='),disp(' '),disp(plim')
  disp(' Blue (outer) set: { v : v = B*u, umin < u < umax }')
  disp('  Feasible virtual control set with constrained allocation')
  disp(' ')
  disp(' Red (inner) set: { v : umin < P*v < umax } where P = pinv(B)')
  disp('  Feasible virtual control set with linear allocation, u = P*v')
  disp(' ')
  disp(sprintf(' Red to blue size ratio: %0.3g%%',ratio*100))
  disp(' ')
  disp('---------------------------------------------------------------')
  
  end