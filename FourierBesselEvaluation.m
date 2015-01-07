function f = FourierBesselEvaluation( nu, c, tol )
% FOURIERBESSELEVALUATION   Evaluate a Fourier-Bessel expansion at (1:N)/N.
% 
% F = FOURIERBESSELEVALUATION( NU, C, TOL ) evaluates the expansion 
% f(r) = sum_{n=1}^N C(n) J_{NU}(j(n) r) at r = 1/N,2/N,...,N/N, where N is 
% length of C, J_{NU}(z) is a Bessel function of parameter NU, and 
% j(n) is the nth positive root of J_0(z). NU should be an integer. 
% The values are return in the vector F with an expected 
% accuracy of ||F||_inf = O(TOL||C||_1), where ||.||_1 is the vector 1-norm
% and ||.||_inf is the absolute maximum vector norm.  
% 
% ALGORITHM: 
%   The algorithm is based on carefully approximating J_{NU}(z) by an 
%   asymptotic expansion for large arguments. The algorithm requires 
%   O( N (log N)^2 / loglog N ) complexity. For more information, see 
%   Section 5 of [1]. 
%
% [1] A. Townsend, A fast analysis-based discrete Hankel transform using 
%     asymptotic expansions, SIAM J. Numer. Anal., submitted, 2015. 
% 
% Author:  
%    Alex Townsend, Dec 14 (originally written)
%    Alex Townsend, Jan 15 (modified)

% The expansion coefficients should be supplied as a column vector: 
[N, m] = size( c ); 
if ( m > 1 ) 
    error('FASTSCHLOMILCHEVALUATION:INPUTS',...
                     'Expansion coefficients should be a column vector.') 
end

% Preallocate output vector: 
f = zeros( N, 1 ); 

% Determine the algorithmic parameters (see Table 4.1 in [1]): 
[M, alpha, beta, P, K, T] = DetermineParameters( nu, N, tol ); 

% We cannot use (5.5) in [1] on the first floor( min(pK, qT) ) columns:
pK = min( exp(1)/(16*pi)*(5.2/tol)^(1/K) + .25, N);            % See (5.6)
qT = min( tol^(-1/(2*T))/(16*pi*factorial(T)^(1/T)) + .25, N); % See (5.7)
c2 = c( 1 : floor( max( pK, qT ) ) ); % Split c into two vectors c1 and c2.
c( 1 : floor( max( pK, qT ) ) ) = 0; 

% We want to compute J_{NU}( r j0^T )*c, where j0 = tildew + b: 
r = (1:N).'./N;   
tildew = ((1:N) -1/4)*pi;   
j0 = besselroots( N );  b = j0 - tildew';

%%%%%%%%%%%%%%%%%%%%%    Bessel Evaluations   %%%%%%%%%%%%%%%%%%%%
% It is cheaper to compute all the Bessel evaluation now, instead of during
% each Schlomilch evaluation.

% Partition parameters (see Figure 4.3 in [1]): 
qo = [1  floor(alpha*beta.^(P:-1:0)*sqrt(N))]; % qo = "Q-odd" for Q_{2p+1}
qe = [floor(alpha*beta.^(0:-1:-P)*sqrt(N)) N]; % qe = "Q-even" for Q_{2p}

% Compute all the Bessel evaluations in J_{NU}^{EVAL} that we will need 
% for direct summation part:  B1 = cell array for horizontal rectangles in
% Figure 4.3, and B2 = cell array for vertical rectangles.
z1 = r(qo(1):qo(2)) * tildew; z2 = r(qe(1)+1:N)*tildew(1:qo(end-P));
B1{1}{1} = besselj(0,z1); B1{1}{2} = besselj(1,z1);
B2{P+1}{1} = besselj(0,z2); B2{P+1}{2} = besselj(1,z2);
for j = 2:nu+2*T+K-3 - max(floor(2*T-2),0) + 1;
    B1{1}{j+1} = 2*(j-1)*B1{1}{j}./z1 - B1{1}{j-1};
    B2{P+1}{j+1} = 2*(j-1)*B2{P+1}{j}./z2 - B2{P+1}{j-1};
end

for p = 2:P+1
    Q2kp1 = qo(p)+1:qo(p+1); cQ2k = 1:qe(end-p+1);      % partition indices
    Q2k = qe(1)+1:qe(p); cQ2kp1 = qo(end-p+1)+1:qo(end-p+2);
    z1 = r(Q2kp1) * tildew(cQ2k); z2 = r(Q2k) * tildew(cQ2kp1);
    B1{p}{1} = besselj(0,z1); B1{p}{2} = besselj(1,z1);
    B2{p-1}{1} = besselj(0,z2); B2{p-1}{2} = besselj(1,z2);
    for j = 2:nu+2*T+K-3 - max(floor(2*T-2),0) + 1;
        B1{p}{j+1} = 2*(j-1)*B1{p}{j}./z1 - B1{p}{j-1};
        B2{p-1}{j+1} = 2*(j-1)*B2{p-1}{j}./z2 - B2{p-1}{j-1};
    end
end

% Now, compute the sum in Section 5.1 (computation has been rearranged for
% greater efficiency): 
for u = 0:2*T+K-3
    
    tt = max(ceil((u-K+1)/2),0):min(floor(u/2),T-1);% inner sum indices
    ss = u-2*tt;                            % Change of variables indicies.
    
    % For Bessel functions with large parameters, we hardly need any terms
    % in the Taylor series. Take account of this: 
    minss = min( ss ); 
    if ( minss > 0 )
       TT = 1:10; 
       tmpT = find( ( eps.^(-1./(2*TT+minss))/(16*pi)./factorial(TT).^(1./(TT+minss)) + .25 ) ...
                                                   < 30, 1, 'first' );
       tt = max(ceil((u-K+1)/2),0):min(floor(u/2),tmpT-1);
       ss = u - 2*tt; 
    end
    if ( isempty( tt ) ) 
        continue
    end
    
    
    Dr = r.^u; Db = ( b /2 ).^u;     % Diagonal matrices in (5.5) 
    
    % Compute inner sum (over t) with a single Schlomilch evaluation: 
    u1 = unique(abs(nu-ss)); u2 = unique(nu+ss); % Which NU's do we need?
    NNU = unique([u1(:) ; u2(:)]);
    tmp3 = VectorizedSchlomilch( NNU, Db.*c, -1/4, r, tildew, ...
                                                M, alpha, beta, P, B1, B2);
    
    % Now figure out how this relates to J_{NU-u+2t} and J_{NU+u-2t}: 
    idx = (repmat( NNU, size(ss) ) == repmat( abs(nu-ss), size(NNU) ));
    [ii, ignored] = ind2sub(size(NNU).*size(ss),find(idx));
    tmp1 = tmp3(:,ii);
    idx = (repmat( NNU, size(ss) ) == repmat( nu+ss, size(NNU) ));
    [ii, ignored] = ind2sub(size(NNU).*size(ss),find(idx));
    tmp2 = tmp3(:,ii);

    % Now compute the u term in rearranged computation: 
    s1 = (-1).^((nu-ss).*((nu-ss)<0)); s2 = (-1).^ss;
    const = (-1).^tt./factorial(tt)./gamma(u-tt+1);
    if ( mod(u,2) == 0 )
        const(end) = const(end)/2;    % Double prime on the sum.
    end 
    tmp = (tmp1 * diag( s1 ) + tmp2 * diag( s2 )) * diag( const );
    f = f + Dr .* sum(tmp,2);                   % Do inner sum.
    
end

% Compute f2 = J_{NU}( rj0^T )*c2 using direct summation and add to f1:
f = f + besselj( nu, r * j0( 1 : floor( max( pK, qT ) ) )' ) * c2;

end

function f = VectorizedSchlomilch( nu, c, shift, r, w, ...
                                                M, alpha, beta, P, B1, B2 )
% VECTORIZEDSCHLOMILCH   Vectorized evaluation of Schlomilch expansions
% 
% F = VECTORIZEDSCHLOMILCH( NU, C, SHIFT, <PARAMETERS> ) evaluates the 
% expansion f(r) = sum_{n=1}^N C(n) J_{NU}((n+SHIFT) pi r) at r = 1/N,2/N,...,N/N, 
% where N is length of C and J_{NU}(z) is a Bessel function of parameter NU. 
% NU can be vector of integers. 
%
% [1] A. Townsend, A fast analysis-based discrete Hankel transform using 
%     asymptotic expansions, SIAM J. Numer. Anal., submitted, 2015. 
% 
% Author:  
%    Alex Townsend, Dec 14 (originally written)
%    Alex Townsend, Jan 15 (modified)

% Length of Schlomilch expansion
N = numel( c ); NUK = numel( nu );

% Preallocate output vector: 
f = zeros( N, NUK ); 

% Partition parameters (see Figure 4.3 in [1]): 
qo = [1  floor(alpha*beta.^(P:-1:0)*sqrt(N))]; % qo = "Q-odd" for Q_{2p+1}
qe = [floor(alpha*beta.^(0:-1:-P)*sqrt(N)) N]; % qe = "Q-even" for Q_{2p}

%%%%%%%%%  Calculate the J_{NU}^{EVAL}( rw^T )*c in (4.9) in [1] %%%%%%%%
% See Figure 4.3: 
for j = 1:NUK    % top horiz partition and far-left vert partition
    f(qo(1):qo(2),j) = f(qo(1):qo(2),j) + B1{1}{nu(j)+1} * c; 
    f(qe(1)+1:N,j) = f(qe(1)+1:N,j) + B2{P+1}{nu(j)+1} * c(1:qo(2));
end

for p = 2:P+1
    Q2kp1 = qo(p)+1:qo(p+1); cQ2k = 1:qe(end-p+1);      % partition indices
    Q2k = qe(1)+1:qe(p); cQ2kp1 = qo(end-p+1)+1:qo(end-p+2);
    for j = 1:NUK
        f(Q2kp1,j) = B1{p}{nu(j)+1} * c(cQ2k);
        f(Q2k,j) = f(Q2k,j) + B2{p-1}{nu(j)+1} * c(cQ2kp1);
    end
end

%%%%%%%%%  Calculate the J_{NU}^{ASY}( rw^T )*c in (4.9) in [1] %%%%%%%%

Pc = c; Pc(1:qe(1)) = 0;                    % Q1 * c 
f_asy = VectorizedQJASYQ(nu,Pc,r,w,M,shift);% J_{NU}^{ASY}( rw^T )*(Q1*c)
f_asy(1:qo(end),:) = 0; f = f + f_asy;      % Q1*J_{NU}^{ASY}( rw^T )*(Q1*c)

for p = 1:P
    Pc = c; Pc(1:qe(end-p)) = 0;             % Q2p * c
    f_asy = VectorizedQJASYQ(nu,Pc,r,w,M,shift);% J_{NU}^{ASY}( rw^T )*(Q2p*c)
    f_asy([1:qo(p+1) qo(p+2)+1:end],:) = 0; % Q2p1*J_{NU}^{ASY}( rw^T )*(Q2p*c)
    f = f + f_asy;
       
    Pc = c; Pc([1:qo(end-p) qo(end-p+1)+1:N]) = 0;% Q2p1 * c
    f_asy = VectorizedQJASYQ(nu,Pc,r,w,M,shift);% J_{NU}^{ASY}( rw^T )*(Q2p1*c)
    f_asy(1:qe(p+1),:) = 0;                    % Q2p*J_{NU}^{ASY}( rw^T )*(Q2p1*c)
    f = f + f_asy;
end

end

function F = VectorizedQJASYQ( nu, c, r, w, M, shift )
% VECTORIZEDQJASYQ      Compute J_{NU}( rw^T ) * c using asymptotic expansion. 
% 
% F = VECTORIZEDQJASYQ( NU, C, R, W, M, SHIFT) computes the matrix-vector 
% product J_{NU}( R*W^T ) * C using 2M terms of the asymptotic expansion for 
% Bessel functions of large arguments. NU can be a vector. 
% 
% This code is based on (4.5) in [1]. 
% 
% Author:  
%    Alex Townsend, Dec 14 (originally written)
%    Alex Townsend, Jan 15 (modified)
 
w = w';     % Make w a column vector 
NUK = numel( nu ); 

% Preallocate output: 
F = zeros( numel(r), NUK ); 

for j = 1:NUK
    % Constants d1 and d2 in (4.3). Also see Remark 4.1: 
    d1(:,j) = cos( -r*shift*pi + (2*nu(j) + 1)*pi/4);  
    d2(:,j) = sin( -r*shift*pi + (2*nu(j) + 1)*pi/4); 
end

ak = ones(1, NUK);    % a_0(nu) = 1. 
rk = 1./sqrt(r);      % 1/sqrt(z)  term
wk = c./sqrt(w);      % 1/sqrt(z)  term

% First term of ASY for each NU: 
QDCT1Q = Qdct1Q( wk ); QDST1Q = Qdst1Q( wk );
for j = 1:NUK
    F(:,j) = ( ( QDCT1Q.*d1(:,j) + QDST1Q.*d2(:,j)).*rk );
end

sgn = -1;
for k = 1:2:2*M
    % Terms containing sin(mu): 
    for j = 1:NUK
        ak(j) = (4*nu(j)^2 - (2*k-1)^2) * ak(j) / k / 8;   % a_{2m+1}(nu)
    end
    rk = rk./r; wk = wk./w;
    
    QDCT1Q = Qdct1Q( wk ); QDST1Q = Qdst1Q( wk );
    for j = 1:NUK
        F(:,j) = F(:,j) + sgn * ak(j) * ((-QDCT1Q.*d2(:,j) +...
                                                   QDST1Q.*d1(:,j) ).*rk);
    end
    
    % Terms containing cos(mu):
    if ( k < 2*M-2 )
        for j = 1:NUK
            ak(j) = (4*nu(j)^2 - (2*k+1)^2) * ak(j) / (k+1) / 8; % a_{2m}(nu)
        end
        rk = rk./r; wk = wk./w;
        
        QDCT1Q = Qdct1Q( wk ); QDST1Q = Qdst1Q( wk );
        for j = 1:NUK
            F(:,j) = F(:,j) + sgn * ak(j) * (( QDCT1Q.*d1(:,j) +...
                                                   QDST1Q.*d2(:,j) ).*rk);
        end
    end
    sgn = -sgn;
end

F = sqrt(2/pi) * F;           % constant out the front

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
