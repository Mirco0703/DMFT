function anderson_green_inverted=andGreen(iw_n,x,mu)
% inverted anderson greensfunction G_0^And^-1(iw_n) on matsubara frequencies 
%
%   Args:
%       iw_n:                    matsubara frequencies
%       x:                       v_l and eps_l, parameters for the anderson greensfunction
%       mu:                      chemical potential
%
%   Returns:
%       anderson_green_inverted: inverted anderson greensfunction

anderson_green_inverted=zeros(1,length(iw_n));
[v_l,eps_l]=para(x);            % gaining full parameters out of short x vector
for ii=1:length(iw_n)
    anderson_green_inverted(ii)=iw_n(ii)+mu-sum(v_l.^2./(iw_n(ii)-eps_l));
end