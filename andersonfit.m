function xmin = andersonfit(mu, n_s, iw_n, selfen, eps_k)
% Fitting of an anderson green function to the bath green function using
% quasi-newton-algorithm of built-in fminunc and the in chi.m defined cost 
% function
%
%   Args:
%       mu:     chemical potential
%       ns:     number of anderson sites
%       w:      matsubara frequencies iw_n
%       selfen: local selfenergy Sigma(iw_n)
%       eps:    energy eigenvalues epsilon(k)
%
%   Returns:
%       xmin:   optimized parameters for the Anderson greensfunction, use
%               para.m to get v_l and epsilon_l from them

fprintf('\n->Fitting of anderson model<-\n')
bath_green_inverted=bathGreen(iw_n,mu,eps_k,selfen);

%define cost function as function of only x for fminunc 
f=@(x) chi(x,iw_n,mu,bath_green_inverted);

%starting values with small randomness
if rem(n_s,2)
	x0=[0.1.*ones(1,ceil(n_s/2)) 0.1+0.1.*rand(1,floor(n_s/2))];
else
	x0=[0.1.*ones(1,n_s/2) 0.1+0.1.*rand(1,n_s/2)];
end

%optimization and printing
op=optimoptions('fminunc','Algorithm','quasi-newton','Display','off');
[xmin,cost]=fminunc(f,x0,op);
[v_l,eps_l]=para(xmin);
fprintf('\n v_l:         eps_l:\n')
fprintf('%10f %10f\n', [v_l; eps_l])
fprintf('\n Chi²: %10f\n',cost)
