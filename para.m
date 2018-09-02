function [vl, epsl] = para(x)
% returns parameters of the anderson impurity model using symmetry (v_l
% in pairs the same, epsilon_l with different sign), for example for ns=4:
% x1 -> v_1 v_1
% x2 -> v_2 v_2
% x3 -> eps_1 -eps_1
% x4 -> eps_2 -eps_2
% For odd number of elements in x, the middle is the last v_l,
% corresponding to epsilon_l=0 for odd number of bath sites
%
%   Args:
%       x:  independent parameters, in first half v_l, in second epsilon_l
%
%   Returns:
%       vl:     v_l
%       epsl:   epsilon_l

ns=length(x);
vl=zeros(1,ns);
epsl=zeros(1,ns);
for ii=1:floor(ns/2)
    vl(2*ii-1)=x(ii);
    vl(2*ii)=x(ii);
    epsl(2*ii-1)=x(ii+ceil(ns/2));
    epsl(2*ii)=-x(ii+ceil(ns/2));
end
if rem(ns,2)
	vl(ns)=x(ceil(ns/2));
end
