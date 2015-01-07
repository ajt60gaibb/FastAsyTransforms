function f = directSummation( nu, r, w, c ) 
% DIRECTSUMMATION    Compute J_{NU}( rw^T )*c using direct summation
% 
%  F = DIRECTSUMMATION( NU, R, W, C ) computes the matrix-vector produce 
%  J_{NU}( rw^T )*c by naively evaluating Bessel functions and performing
%  the matrix-vector product in O(N^2) cost. 
%  
% Author:  
%    Alex Townsend, Dec 14 (originally written)
%    Alex Townsend, Jan 15 (modified)

% Do not form a matrix since for N>5000 I have 
% memory issues: 
f = zeros(numel(r)); 
for k = 1:numel(r)
    f(k) = besselj( nu, r(k)*w )*c;
end

end 