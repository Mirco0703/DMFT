function local_green=locGreen(iw_n,mu,eps_k,selfen)
% calculate the local Green function on matsubara frequencies
%
%   Args:
%       iw_n:        Matsubara frequencies
%       mu:          chemical potential
%       eps_k:       dispersion in tight-binding
%       selfen:      self-energy on Matsubara frequenies
%
%   Returns:
%       local_green: local Green function

local_green=zeros(1,length(iw_n));
for ii=1:length(iw_n)
    % sum(sum()) is the sum over k, because eps_k is the only non-scalar
    % double sum because eps_k is a matrix
    local_green(ii)=sum(sum(1./(iw_n(ii)-eps_k+mu-selfen(ii))))./numel(eps_k);
end