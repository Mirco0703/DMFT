function epsk = tightbinding(a,t)
% Calculating dispersion and density of states for quadratic lattice in
% tight-binding-approximation in the first brillouin zone
%
%   Args:
%       a:  lattice constant
%       t:  hopping parameter
%
%   Returns:
%       epsk:    dispersion

% lattice vectors
v1=[a; 0];
v2=[0; a];

% size of first brillouin zone
g=2*pi/a;

% define k vectors in first BZ
grid_size=1000;
k=zeros(grid_size,grid_size,2);
for ii=1:grid_size
    for jj=1:grid_size
        k(ii,jj,1)=(ii)/grid_size*g;
        k(ii,jj,2)=(jj)/grid_size*g;
    end
end

% calculating dispersion
epsk=zeros(grid_size);
for ii=1:grid_size
    for jj=1:grid_size
        epsk(ii,jj)=t*(exp(1i*squeeze(k(ii,jj,:))'*v1) + exp(-1i*squeeze(k(ii,jj,:))'*v1) + exp(1i*squeeze(k(ii,jj,:))'*v2) + exp(-1i*squeeze(k(ii,jj,:))'*v2));
    end
end



