function f = DHT( c, tol )
% DHT   Compute the discrete Hankel transform
% 
% F = DHT( C, TOL ) computes the sums sum_{n=1}^N C(n)J_{0}(j(n)j(k)/j(N)), 
% where N is length of C, J_{0}(z) is a Bessel function, and 
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

% Whoops: Originally coded up with the wrong size. Minor fix: 
c = [ c ; 0 ]; N = N + 1; 

% Preallocate output vector: 
f = zeros( N, 1 ); 

% Determine the algorithmic parameters (see Table 4.1 in [1]): 
[M, alpha, beta, P, K, T] = DetermineParameters( 0, N, tol ); 

% We cannot use (5.5) in [1] on the first 1.01floor( min(pK, qT) ) rows 
% and columns:
pK = exp(1)/(16*pi)*(5.2/tol)^(1/K) + .25;            % See (5.6)
qT = tol^(-1/(2*T))/(16*pi*factorial(T)^(1/T)) + .25; % See (5.7)
max_cut = min( floor( 1.01*max( pK, qT ) ), N);
c2 = c( 1 : max_cut );             % Split c into two vectors c1 and c2...
Pc = c; Pc( 1 : max_cut ) = 0;     % and zero out first few entries.

% We want to compute J_{0}( j0 j0^T/j0(N+1) )*c, where 
% j0/j0(N+1) = tilde_r + eN and j0 = tildew + b: 
j0 = besselroots( N ); ratios = j0 ./ j0( N );
tilde_r = (4*(1:N)'-1)./(4*N-1);
eN = ratios - tilde_r;         %  j0/j0(N+1) = tilde_r + eN
tildew = ((1:N) -.25)*pi;
b = j0 - tildew';              %  j0 = tildew + b

%%%%%%%%%%%%%%%%%%%%%    Bessel Evaluations   %%%%%%%%%%%%%%%%%%%%
% It is cheaper to compute all the Bessel evaluation now, instead of during
% each Schlomilch evaluation.

% Partition parameters (see Figure 4.3 in [1]): 
qo = [1  floor(alpha*beta.^(P:-1:0)*sqrt(N))]; % qo = "Q-odd" for Q_{2p+1}
qe = [floor(alpha*beta.^(0:-1:-P)*sqrt(N)) N]; % qe = "Q-even" for Q_{2p}

% Compute all the Bessel evaluations in J_{NU}^{EVAL} that we will need 
% for direct summation part:  B1 = cell array for horizontal rectangles in
% Figure 4.3, and B2 = cell array for vertical rectangles.
z1 = tilde_r(qo(1):qo(2)) * tildew; 
z2 = tilde_r(qe(1)+1:N)*tildew(1:qo(end-P));
B1{1}{1} = besselj(0,z1); B1{1}{2} = besselj(1,z1);
B2{P+1}{1} = besselj(0,z2); B2{P+1}{2} = besselj(1,z2);
for j = 2:2*T+2*K-4 - max(floor(2*T-2),0) + 1;
    B1{1}{j+1} = 2*(j-1)*B1{1}{j}./z1 - B1{1}{j-1};
    B2{P+1}{j+1} = 2*(j-1)*B2{P+1}{j}./z2 - B2{P+1}{j-1};
end

for p = 2:P+1
    Q2kp1 = qo(p)+1:qo(p+1); cQ2k = 1:qe(end-p+1);      % partition indices
    Q2k = qe(1)+1:qe(p); cQ2kp1 = qo(end-p+1)+1:qo(end-p+2);
    z1 = tilde_r(Q2kp1) * tildew(cQ2k); z2 = tilde_r(Q2k) * tildew(cQ2kp1);
    B1{p}{1} = besselj(0,z1); B1{p}{2} = besselj(1,z1);
    B2{p-1}{1} = besselj(0,z2); B2{p-1}{2} = besselj(1,z2);
    for j = 2:2*T+2*K-4 - max(floor(2*T-2),0) + 1;
        B1{p}{j+1} = 2*(j-1)*B1{p}{j}./z1 - B1{p}{j-1};
        B2{p-1}{j+1} = 2*(j-1)*B2{p-1}{j}./z2 - B2{p-1}{j-1};
    end
end

% Now, compute the sum in (6.2) (computation has been rearranged for
% greater efficiency): 
for u = 0:2*T+K-3
    
    tt = max(ceil((u-K+1)/2),0):min(floor(u/2),T-1);% inner sum indices
    ss = u-2*tt;                           % Change of variables indicies.
    
    % For Bessel functions with large parameters, we hardly need any terms
    % in the Taylor series. Take account of this: 
    minss = min( ss ); 
    if ( minss > 0 )
       TT = 1:10; 
       tmpT = find( tol.^(-1./(2*TT+minss))/(16*pi)./factorial(TT).^(1./(TT+minss) + .25 ) ...
                                                   < 30, 1, 'first' );
       tt = max(ceil((u-K+1)/2),0):min(floor(u/2),tmpT-1);
       ss = u - 2*tt; 
    end
    if ( isempty( tt ) ) 
        continue
    end
    
    DeN = (eN*(N-1/4)).^u;                 % Diagonal matrices in (6.2)
    Dj0 = (j0/2/(N-1/4)).^u; 
    
    % Now compute the u term in rearranged computation: 
    CONST = 2*(-1).^tt./factorial(tt)./gamma(ss+tt+1);
    if ( mod( u, 2 )==0 ), CONST(end) = CONST(end)/2; end
    SGN = (-1).^ss;
    tmp = VectorizedFourierBessel( ss, (Dj0.*Pc), tilde_r, tildew, b, B1, B2, ...
                                                  M, alpha, beta, P, K, T);
    f = f + bsxfun(@times, tmp * ( SGN.*CONST )', DeN); % Do inner sum.
end

% Compute f2 = J_{NU}( j0j0^T/j0(N+1) )*c2 using direct summation and 
% add to f1:
f = f + besselj( 0, ratios*j0( 1: max_cut )' ) * c2;
% Compute the columns of f that we could not use (6.2) on: 
f(1:max_cut) = besselj(0, ratios(1:max_cut)*j0') * c;

% Whoops: Originally coded up with the wrong size. Minor fix: 
f(N) = []; 

end



function f = VectorizedFourierBessel( nu, c, tilde_r, tildew, b, B1, B2, ...
                                                   M, alpha, beta, P, K, T)
% VECTORIZEDFOURIERBESSEL   Evaluate a Fourier-Bessel expansion at (1:N)/N.
% 
% F = VECTORIZEDFOURIERBESSEL( NU, C, TOL ) evaluates the expansion 
% f(r) = sum_{n=1}^N C(n) J_{NU}(j(n) r) at 
% r = 3/(4N+3),7/(4N+3),...,(4N-1)/(4N+3), where N is 
% length of C. NU can be a vector of integers. 
%
% [1] A. Townsend, A fast analysis-based discrete Hankel transform using 
%     asymptotic expansions, SIAM J. Numer. Anal., submitted, 2015. 
% 
% Author:  
%    Alex Townsend, Dec 14 (originally written)
%    Alex Townsend, Jan 15 (modified)

%%%%%%%%%%%%%%% Calculate J^{ASY}P*c = J^{ASY}*Pc %%%%%%%%%%%%%%%%
f = zeros( length(tilde_r), numel(nu) ); % f_asy = contribution to output from P^TJ^{ASY}P
for u = 0:2*T+K-3
    
    tt = max(ceil((u-K+1)/2),0):min(floor(u/2),T-1);% inner sum indices
    ss = u-2*tt;                           % Change of variables indicies.
    
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
    
    Dr = tilde_r.^(u); Db = (b/2).^(u);     % Diagonal matrices in (5.5) 

    % Compute inner sum (over t) with a single Schlomilch evaluation:
    NU = repmat( nu', size(ss) );   
    SS = repmat( ss, size(nu') );
    u1 = unique(abs(NU-SS)); u2 = unique(NU+SS);  % Which NU's do we need?
    NNU = unique([u1(:) ; u2(:)]);
    tmp3 = VectorizedSchlomilchEvaluation( NNU, Db.*c, tilde_r, tildew,...
                                                M, alpha, beta, P, B1, B2);
    
    % For each NU, figure out how this relates to J_{NU-u+2t} and J_{NU+u-2t}:
    for j = 1:numel(nu)  
        idx = (repmat( NNU, size(ss) ) == repmat( abs(nu(j)-ss), size(NNU) ));
        [ii, ignored] = ind2sub(size(NNU).*size(ss),find(idx));
        tmp1 = tmp3(:,ii);
        idx = (repmat( NNU, size(ss) ) == repmat( nu(j)+ss, size(NNU) ));
        [ii, ignored] = ind2sub(size(NNU).*size(ss),find(idx));
        tmp2 = tmp3(:,ii);
        
        %For each NU, compute the u term in rearranged computation: 
        s1 = (-1).^((nu(j)-ss).*((nu(j)-ss)<0)); sgn = (-1).^ss;
        const = (-1).^tt./gamma(tt+1)./gamma(u-tt+1);
        if ( mod(u,2) == 0), const(end) = const(end)/2; end
        tmp = (tmp1 * diag( s1 ) + tmp2 * diag( sgn )) * const'; 
        f(:,j) = f(:,j) + Dr .* tmp;              % Do inner sum.
    end
end

end


function f = VectorizedSchlomilchEvaluation( nu, c, tilde_r, tildew,...
                                                M, alpha, beta, P, B1, B2 )
% VECTORIZEDSCHLOMILCH   Vectorized evaluation of Schlomilch expansions
% 
% F = VECTORIZEDSCHLOMILCH( NU, C, SHIFT, <PARAMETERS> ) evaluates the 
% expansion f(r) = sum_{n=1}^N C(n) J_{NU}((n+SHIFT) pi r) at r = 1/N,2/N,...,N/N, 
% where N is length of C and J_{NU}(z) is a Bessel function of parameter NU. 
% NU can be a vector of integers. 
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
f_asy = VectorizedQJASYQ(nu,Pc,tilde_r,tildew,M);% J_{NU}^{ASY}( rw^T )*(Q1*c)
f_asy(1:qo(end),:) = 0; f = f + f_asy;      % Q1*J_{NU}^{ASY}( rw^T )*(Q1*c)

for p = 1:P
    Pc = c; Pc(1:qe(end-p)) = 0;             % Q2p * c
    f_asy = VectorizedQJASYQ(nu,Pc,tilde_r,tildew,M);% J_{NU}^{ASY}( rw^T )*(Q2p*c)
    f_asy([1:qo(p+1) qo(p+2)+1:end],:) = 0; % Q2p1*J_{NU}^{ASY}( rw^T )*(Q2p*c)
    f = f + f_asy;
       
    Pc = c; Pc([1:qo(end-p) qo(end-p+1)+1:N]) = 0;% Q2p1 * c
    f_asy = VectorizedQJASYQ(nu,Pc,tilde_r,tildew,M);% J_{NU}^{ASY}( rw^T )*(Q2p1*c)
    f_asy(1:qe(p+1),:) = 0;                    % Q2p*J_{NU}^{ASY}( rw^T )*(Q2p1*c)
    f = f + f_asy;
end

end

function F = VectorizedQJASYQ( nu, c, r, w, M )
% VECTORIZEDQJASYQ    Compute J_{NU}( rw^T ) * c using asymptotic expansion. 
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

w = w'; shift = (w(1)/pi - 1); 
N = numel(r); NUK = numel( nu ); 

% Preallocate output: 
F = zeros( N, NUK ); 

d1 = zeros( N, NUK); d2 = d1; 
for j = 1:NUK
    % Constants d1 and d2 in (4.3). Also see Remark 4.1: 
    d1(:,j) = cos( -r*shift*pi + nu(j)/2*pi + pi/4);
    d2(:,j) = sin( -r*shift*pi + nu(j)/2*pi + pi/4);
end

ak = ones(1, NUK);     % a_0(nu) = 1. 
rk = 1./sqrt(r);       % 1/sqrt(z)  term in (3.1)
wk = c./sqrt(w);       % 1/sqrt(z)  term in (3.1)

% twiddle factor to speed up DCTs and DSTs
t1 = exp(-1i*(1:N)'*pi/(4*N-1)); wk = t1.*wk; 

% First term of ASY for each NU: 
[QDCT1Q, QDST1Q] = DCTs( wk ); 
for j = 1:NUK
    F(:,j) = ( ( QDCT1Q.*d1(:,j) + QDST1Q.*d2(:,j)).*rk );
end

sgn = -1;
for k = 1:2:2*M
    % Terms containing sin(mu) in (3.1): 
    rk = rk./r; wk = wk./w;
    for j = 1:NUK
        ak(j) = (4*nu(j)^2 - (2*k-1)^2) * ak(j) / k / 8;
    end
    [QDCT1Q, QDST1Q] = DCTs( wk ); 
    for j = 1:NUK
        F(:,j) = F(:,j) + sgn * ak(j) * ((-QDCT1Q.*d2(:,j) + QDST1Q.*d1(:,j) ).*rk);
    end
    
    % Terms containing cos(mu) in (3.1):
    if ( k < 2*M-2 )
        for j = 1:NUK
            ak(j) = (4*nu(j)^2 - (2*k+1)^2) * ak(j) / (k+1) / 8;
        end
        rk = rk./r; wk = wk./w;
        [QDCT1Q, QDST1Q] = DCTs( wk ); 
        for j = 1:NUK
            F(:,j) = F(:,j) + sgn * ak(j) * (( QDCT1Q.*d1(:,j) + QDST1Q.*d2(:,j) ).*rk);
        end
    end
    sgn = -sgn;
end

F = sqrt(2/pi) * F;               % constant out the front

end

function [Y1, Y2] = DCTs( u ) 
% DCTs    compute the DCT and DST of type 1 together
% 
%  Y = DCTs( c ) compute the matrix-vector product Q C_{N+1}^I Q * c,
%  where C_{N+1}^I is the DCT matrix of size (N+1)x(N+1) and Q = [0 | I_N]
%  and [S_{N-1}^I * c(1:end-1) ; 0], where S_{N-1}^I is the DST matrix 
%  of size (N-1)x(N-1).
% 
% Author:  
%    Alex Townsend, Dec 14 (originally written)
%    Alex Townsend, Jan 15 (modified)

% Highly optimized code: 
N = numel( u ); 
Y2 = fft(vertcat(0,u,zeros( 3*N - 2, 1))); 
Y1 = real(Y2(2:2:2*N));
Y2 = -imag(Y2(2:2:2*N));

end