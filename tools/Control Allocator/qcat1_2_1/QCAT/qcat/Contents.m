% QCAT - Quadratic Programming Control Allocation Toolbox
% Version 1.2.1 25-Aug-2004
%
% Author:    Ola Härkegård, ola@isy.liu.se
% Download:  http://www.control.isy.liu.se/~ola/qcat
% -------------------------------------------------------
% Information.
%   qcatdoc       - Open QCAT documentation in browser.
%
% QP based control allocation
%   sls_alloc     - Active set, sequential least squares.
%   wls_alloc     - Active set, weighted least squares.
%   wlsc_alloc    - C implementation of wls_alloc.
%   mls_alloc     - Active set, minimal least squares.
%   ip_alloc      - Interior point method.
%   cgi_alloc     - Cascaded generalized inverses method.
%   fxp_alloc     - Fixed-point method.
% 
% Direct allocation.
%   dir_alloc     - Direct control allocation.
% 
% Dynamic control allocation.
%   dca           - Design dynamic control allocation filter.
%
% Simulation.
%   qp_sim        - Response of static QP allocator.
%   dir_sim       - Response of direct allocator.
%   dyn_sim       - Response of dynamic allocator.
%   qcatlib       - Simulink blocks for QP and dynamic allocation.
%
% Constrained allocation vs linear control.
%   lq2ca         - Extract allocator from LQ controller.
%   ca2lq         - Merge allocator with LQ controller.
%   vview         - View feasible virtual control set.
%   vview_demo    - Demo of vview.
%
% Other.
%   iscoplanar    - Test for coplanar controls.
%   wpinv         - Weighted pseudoinverse.
%
% Data files.
%   admiredata    - Admire test data.
%   f18data       - F-18 based test data.