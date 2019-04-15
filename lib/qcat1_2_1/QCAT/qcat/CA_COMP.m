function OUT = CA_COMP(v,u_ant,W,us,B,plim,rlim,T,Wv,W1,W2,alg, gammma,imax,k,m,tol,CA_type)
%% Completed Control Allocator, using all algorithms of CA (qp_ca_COMP) of 
%  normal and dynamic CA (dyn_ca_COMP)
%
% Author: David Torres
% Based on work of Ola Härkegård
% Date: 05/06/2014
% 
% See also: dyn_ca_COMP, qp_ca_COMP QCATDOC ,  QP_SIM , WLS_ALLOC, WLSC_ALLOC, FXP_ALLOC, SLS_ALLOC, MLS_ALLOC, CGI_ALLOC, IP_ALLOC .

switch lower(CA_type)
    case 'dyn'
        OUT = dyn_ca_COMP(v,u_ant,W,us,B,plim,rlim,T,Wv,W1,W2,alg, gammma,imax,k,m,tol);
    case 'qp'
        OUT = qp_ca_COMP(v,u_ant,W,us,B,plim,rlim,T,Wv,W1,alg, gammma,imax,k,m,tol);
end


end

