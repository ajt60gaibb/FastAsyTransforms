function f = SlowSchlomilchEvaluation( nu, c, tol, shift )
% SLOWSCHLOMILCHEVALUATION    Evaluate a Schlomilch expansion at (1:N)/N.
% 
% F = SLOWSCHLOMILCHEVALUATION( NU, C, TOL ) evaluates the expansion 
% f(r) = sum_{n=1}^N C(n) J_{NU}(n pi r) at r = 1/N,2/N,...,N/N, where N is 
% length of C and J_{NU}(z) is a Bessel function of parameter NU. NU should 
% be an integer. The values are return in the vector F with an expected 
% accuracy of ||F||_inf = O(TOL||C||_1), where ||.||_1 is the vector 1-norm
% and ||.||_inf is the absolute maximum vector norm. 
% 
% F = SLOWSCHLOMILCHEVALUATION( NU, C, TOL, SHIFT ) evaluates the expansion 
% f(r) = sum_{n=1}^N C(n) J_{NU}((n+SHIFT) pi r) at r = 1/N,2/N,...,N/N.  
% 
% ALGORITHM: 
%   The algorithm is based on carefully approximating J_{NU}(z) by an 
%   asymptotic expansion for large arguments. The algorithm requires 
%   O( N^{3/2} ) complexity. For more information, see Section 4.1 of [1]. 
%
% [1] A. Townsend, A fast analysis-based discrete Hankel transform using 
%     asymptotic expansions, SIAM J. Numer. Anal., submitted, 2015. 
% 
% Author:  
%    Alex Townsend, Dec 14 (originally written)
%    Alex Townsend, Jan 15 (modified)

% SHIFT defaults to 0: 
if ( nargin < 4 )
    shift = 0; 
end

% The expansion coefficients should be supplied as a column vector: 
[N, m] = size( c ); 
if ( m > 1 ) 
    error('FASTSCHLOMILCHEVALUATION:INPUTS',...
                     'Expansion coefficients should be a column vector.') 
end

% Determine the algorithmic parameters (see Table 4.1 in [1]): 
[M, alpha] = DetermineParameters( nu, N, tol );

% We want to compute J_{NU}( rw^T )*c, where
r = (1:N)' ./ N; 
w = ( (1:N) + shift ) * pi;

% Partition parameter (see Figure 4.2 in [1]): 
po = ceil( alpha * sqrt( N ) );

%%%%%%%%%%%%%%%%%%%%%% Calculate J^{EVAL} %%%%%%%%%%%%%%%%%%%%%% 
left = besselj( nu, r(po+1:end)*w(1:po) ) * c(1:po);
top = besselj( nu, r(1:po)*w ) * c;
f_eval = [top ; left];

%%%%%%%%%%%%%%%%%%%% Calculate Q1J^{ASY}Q1*c %%%%%%%%%%%%%%%%%%%%% 
Pc = c; Pc(1:po) = 0;           % Calculate P*c
f_asy = QJASYQ( nu, Pc, r, w, M, shift );
f_asy(1:po) = 0; 

%%%%%%%%%%%%%%% Add the two contributions together %%%%%%%%%%%%%%% 
f = f_asy + f_eval; 

end

function F = QJASYQ( nu, c, r, w, M, shift ) 
% QJASYQ      Compute J_{NU}( rw^T ) * c using asymptotic expansion. 
% 
% F = QJASYQ( NU, C, R, W, M, SHIFT) computes the matrix-vector product 
% J_{NU}( R*W^T ) * C using 2M terms of the asymptotic expansion for 
% Bessel functions of large arguments. 
% 
% This code is based on (4.5) in [1]. 
% 
% Author:  
%    Alex Townsend, Dec 14 (originally written)
%    Alex Townsend, Jan 15 (modified)

w = w';            % Make w a column vector 
nusq = 4*nu^2;    

% Constants d1 and d2 in (4.3). Also see Remark 4.1: 
d1 = cos( -r*shift*pi + (2*nu + 1)*pi/4);  
d2 = sin( -r*shift*pi + (2*nu + 1)*pi/4); 

ak = 1;              % a_0(nu) = 1. 
rk = 1./sqrt(r);     % 1/sqrt(z)  term
wk = c./sqrt(w);     % 1/sqrt(z)  term

J0 = ( ( Qdct1Q( wk ).*d1 + Qdst1Q( wk ).*d2).*rk ); % first term of ASY
sgn = -1;
for m = 1:2:2*M
    % Terms containing sin(mu): 
    ak = (nusq - (2*m-1)^2) * ak / m / 8;        % a_{2m+1}(nu)
    rk = rk./r; wk = wk./w;
    J0 = J0 + sgn * ak * ((-Qdct1Q( wk ).*d2 + Qdst1Q( wk ).*d1 ).*rk);
    
    % Terms containing cos(mu): 
    if ( m < 2*M-2 )
        ak = (nusq - (2*m+1)^2) * ak / (m+1) / 8; % a_{2m}(nu)
        rk = rk./r; wk = wk./w;
        J0 = J0 + sgn * ak * (( Qdct1Q( wk ).*d1 + Qdst1Q( wk ).*d2 ).*rk);  
    end
    sgn = -sgn;
end

F = sqrt(2/pi) * J0;                     % constant out the front

end

function y = Qdct1Q( c ) 
% QDCT1Q    compute the DCT of type 1 
% 
%  Y = QDCT1Q( c ) compute the matrix-vector product Q C_{N+1}^I Q * c,
%  where C_{N+1}^I is the DCT matrix of size (N+1)x(N+1) and Q = [0 | I_N]. 
% 
% Author:  
%    Alex Townsend, Dec 14 (originally written)
%    Alex Townsend, Jan 15 (modified)

    c = [ 0 ; c ]; N = numel( c );
    c( N ) = 2 * c( N ); 
    y = fft( [ c ; c(N-1:-1:2) ]/2 );       % Relate DCT to FFT. 
    y = y( 2:N );
end

function y = Qdst1Q( c ) 
% QDST1Q    compute the DST of type 1 
% 
%  Y = QDST1Q( c ) returns the vector [S_{N-1}^I * c(1:end-1) ; 0],
%  where S_{N-1}^I is the DST matrix of size (N-1)x(N-1). 
% 
% Author:  
%    Alex Townsend, Dec 14 (originally written)
%    Alex Townsend, Jan 15 (modified)

     c = c(1:end-1);
     n = numel( c );
     y = fft( [ 0 ; c ; 0 ; -c(n:-1:1,:)] ); % Relate DST to FFT.
     y = -1i*( y( 2*n+2:-1:n+3, : )/2 );
     y = [y ; 0];  
end