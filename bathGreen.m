function bath_green_inverted = bathGreen(iw_n,mu,eps_k,selfen)
% inverted bath green function on matsubara frequencies
%
%   Args:
%       iw_n:                   Matsubara frequencies
%       mu:                     chemical potential
%       eps_k:                  dispersion in tight-binding
%
%   Returns:
%       bath_green_inverted:    inverted bath green function

local_green=locGreen(iw_n,mu,eps_k,selfen);
bath_green_inverted=(1./local_green+selfen);
