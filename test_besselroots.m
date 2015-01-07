function pass = test_besselroots( )
% Test BESSELROOTS()
%
% Author:  
%    Alex Townsend, Jan 15 (originally written)

NN = [1:20 100 1000 2000];
j = 1; tol = 100*eps;
for N = NN
    pass(j) = norm( besselj(0, besselroots( N ) ), inf ) < tol;
    j = j + 1;
end

if ( all( pass ) )
    pass = all( pass );
end

end