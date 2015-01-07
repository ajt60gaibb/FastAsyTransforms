function j0k = besselroots( N )
% BESSELROOTS     Compute the positive roots of J0
%
%  J0K = BESSELROOTS( N )   computes the first N positive roots of J0.
%  The first 20 are tabulated, and the others are computing using an
%  asymptotic expansion.
%
%  Author:
%    Alex Townsend, Dec 14 (originally written)

j0k = zeros( max(N, 20), 1); 

% The first 20 are tabulated:
j0k(1:20) = [ 2.404825557695773
    5.520078110286311
    8.653727912911012
    11.79153443901428
    14.93091770848779
    18.07106396791092
    21.21163662987926
    24.35247153074930
    27.49347913204025
    30.63460646843198
    33.77582021357357
    36.91709835366404
    40.05842576462824
    43.19979171317673
    46.34118837166181
    49.48260989739782
    52.62405184111500
    55.76551075501998
    58.90698392608094
    62.04846919022717];

if N > 20
    %  Use the asymptotic formula (10.21.19) from DLMF:
    nn = 21:N; w = (nn - 1/4)*pi;
    nusq = 0;
    j0k(nn) = w - (nusq-1)/8./w -4*(nusq-1)*(7*nusq-31)/3./(8*w).^3 - ...
        32*(nusq-1)*(83*nusq^2-982*nusq+3779)/15./(8*w).^5 -...
        64*(nusq-1)*(6949*nusq^3-153855*nusq^2+1585743*nusq-6277237)/105./(8*w).^7;
else
    j0k = j0k( 1:N );
end

end