function f = FastSchlomilchEvaluation( nu, c, tol, shift )
% FASTSCHLOMILCHEVALUATION    Evaluate a Schlomilch expansion at (1:N)/N.
% 
% F = FASTSCHLOMILCHEVALUATION( NU, C, TOL ) evaluates the expansion 
% f(r) = sum_{n=1}^N C(n) J_{NU}(n pi r) at r = 1/N,2/N,...,N/N, where N is 
% length of C and J_{NU}(z) is a Bessel function of parameter NU. NU should 
% be an integer. The values are return in the vector F with an expected 
% accuracy of ||F||_inf = O(TOL||C||_1), where ||.||_1 is the vector 1-norm
% and ||.||_inf is the absolute maximum vector norm. 
% 
% F = FASTSCHLOMILCHEVALUATION( NU, C, TOL, SHIFT ) evaluates the expansion 
% f(r) = sum_{n=1}^N C(n) J_{NU}((n+SHIFT) pi r) at r = 1/N,2/N,...,N/N.  
% 
% ALGORITHM: 
%   The algorithm is based on carefully approximating J_{NU}(z) by an 
%   asymptotic expansion for large arguments. The algorithm requires 
%   O( N (log N)^2 / loglog N ) complexity. For more information, see 
%   Section 4.2 of [1]. 
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

% Preallocate output vector: 
f = zeros( N, 1 ); 

% Determine the algorithmic parameters (see Table 4.1 in [1]): 
[M, alpha, beta, P] = DetermineParameters( nu, N, tol );

% We want to compute J_{NU}( rw^T )*c, where
r = (1:N)' ./ N; 
w = ( (1:N) + shift ) * pi;

% Partition parameters (see Figure 4.3 in [1]): 
qo = [1  floor(alpha*beta.^(P:-1:0)*sqrt(N))]; % qo = "Q-odd" for Q_{2p+1}
qe = [floor(alpha*beta.^(0:-1:-P)*sqrt(N)) N]; % qe = "Q-even" for Q_{2p}

%%%%%%%%%  Calculate the J_{NU}^{EVAL}( rw^T )*c in (4.9) in [1] %%%%%%%%
% See Figure 4.3: 
f(qo(1):qo(2)) = besselj( nu, r(qo(1):qo(2))*w) * c;% top horiz partition
f(qe(1)+1:N) = f(qe(1)+1:N) + ...                   % left vert partition
                 besselj( nu, r(qe(1)+1:N)*w(1:qo(2))) * c(1:qo(2));

for p = 2:P+1
    Q2kp1 = qo(p)+1:qo(p+1); cQ2k = 1:qe(end-p+1);      % partition indices
    Q2k = qe(1)+1:qe(p); cQ2kp1 = qo(end-p+1)+1:qo(end-p+2);
    f(Q2kp1) = besselj( nu, r(Q2kp1)*w(cQ2k)) * c(cQ2k);         % kth horz        
    f(Q2k) = f(Q2k) + besselj( nu, r(Q2k)*w(cQ2kp1)) * c(cQ2kp1);% kth vert
end

%%%%%%%%%  Calculate the J_{NU}^{ASY}( rw^T )*c in (4.9) in [1] %%%%%%%%
Pc = c; Pc(1:qe(1)) = 0;                  % Q1 * c 
f_asy = QJASYQ( nu, Pc, r, w, M, shift );    % J_{NU}^{ASY}( rw^T )*(Q1*c)
f_asy(1:qo(end)) = 0; f = f + f_asy;      % Q1*J_{NU}^{ASY}( rw^T )*(Q1*c)

for p = 1:P 
    Pc = c; Pc(1:qe(end-p)) = 0;          % Q2p * c
    f_asy = QJASYQ( nu, Pc, r, w, M, shift );% J_{NU}^{ASY}( rw^T )*(Q2p*c)
    f_asy([1:qo(p+1) qo(p+2)+1:end]) = 0; % Q2p1*J_{NU}^{ASY}( rw^T )*(Q2p*c)
    f = f + f_asy;
    
    Pc = c; Pc([1:qo(end-p) qo(end-p+1)+1:N]) = 0; % Q2p1 * c
    f_asy = QJASYQ( nu, Pc, r, w, M, shift );% J_{NU}^{ASY}( rw^T )*(Q2p1*c)
    f_asy(1:qe(p+1)) = 0;                 % Q2p*J_{NU}^{ASY}( rw^T )*(Q2p1*c)
    f = f + f_asy;
end

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
