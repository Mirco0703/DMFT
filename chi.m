function chi = chi(x,iw_n,mu,bath_green_inverted)
% cost function for fitting the anderson impurity model to the bath
% greensfunction
%
%   Args:
%       iw_n:                   matsubara frequencies
%       x:                      parameters of the anderson impurity model
%       mu:                     chemical potential
%       bath_green_inverted:    inverted bath greensfunction on matsubara frequencies
%
%   Returns:
%       chi:                    cost function Chi^2

chi=sum(abs(bath_green_inverted-andGreen(iw_n,x,mu)).^2)./(1+length(iw_n));