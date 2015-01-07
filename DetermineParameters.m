function [M, alpha, beta, P, K, T] = DetermineParameters( nu, N, tol )
% DETERMINEPARAMETERS     Determine the algorithmic parameters for [1]
% 
% [M, ALPHA, BETA, P] = DetermineParameters( NU, N, TOL ) returns the
% values of M, ALPHA, BETA, and P as given in Table 4.1 in [1]. NU can be a
% vector.
%
% [M, ALPHA, BETA, P] = DetermineParameters( NU, N, TOL ) returns in
% addition the values of K and T used in Section 5 and 6 of [1]. 
%
% [1] A. Townsend, A fast analysis-based discrete Hankel transform using 
%     asymptotic expansions, SIAM J. Numer. Anal., submitted, 2015. 
% 
% Author:  
%    Alex Townsend, Dec 14 (originally written)
%    Alex Townsend, Jan 15 (modified)

% See Table 4.1 for formulas for M, sM, alpha, beta, and P. 
M = max( floor( 0.3*log(1/tol) ), 3);   % Number of terms in ASY expansion. 

% Calculate s_{M,nu}(epsilon) by a fixed point iteration:
sM = 1;
for nuk = nu
    SQRT2 = @(k) 1/2/sqrt(2):1/sqrt(2):k/sqrt(2)-1/2/sqrt(2);
    aM = @(k) prod( (nuk^2/2-SQRT2(k).^2)./(1:k) );
    sMnu = 1;
    for fixedPoint = 1:4
        sMnu = ( sqrt(2)/sqrt(pi)*(abs(aM(2*M)) + ... 
                                abs(aM(2*M+1))/sMnu) / tol )^(1/(2*M+.5));
    end
    sM = max( sM, sMnu );
end

% Partitioning parameters 
alpha = sqrt( sM / pi );
beta = min(3 / log( N ), .8);     % Don't go over beta = 0.8.

% Number of partitions: 
minc = 30; 
P = max( ceil( log( minc/alpha/sqrt(N) ) / log( beta ) ), 0); 

% Number of terms in Neumann addition formula: 
K = ceil( log(5.2/tol) / log( 16*pi*(minc - .25)/exp(1) ) );

% Number of terms in Taylor series expansion: 
TT = 1:10; 
T = find( ( tol.^(-1./(2*TT))/(16*pi)./factorial(TT).^(1./TT) + .25 ) ...
                                                   < minc, 1, 'first' );

end