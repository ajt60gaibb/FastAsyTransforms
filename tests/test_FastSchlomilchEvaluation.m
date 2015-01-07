function pass = test_FastSchlomilchEvaluation( ) 
% Test FASTSCHLOMILCHEVALUATION()
% 
% Author:  
%    Alex Townsend, Jan 15 (originally written)

NN = [100 113 200 234 1000 1001 1013 2089];  % Lots of N. 
TOL = [1e-15 1e-8 1e-3]; 
NU = 0:5; 
pass = ones( numel(NN), numel(NU), numel(TOL) );
j1 = 1;
for tol = TOL
    j2 = 1; 
    for nu = NU
        j3 = 1; 
        for N = NN 
            
            % Expansion coefficients:
            c = randn(N, 1); 
            
            % Our O( N(log N)^2/loglog N ) algorithm: 
            f = FastSchlomilchEvaluation( nu, c, tol );

            % Direct summation: 
            r = (1:N)'./N; w = (1:N)*pi; 
            exact = besselj( nu, r*w )*c;
            
            % Compare error: 
            if ( norm( exact - f, inf ) > tol*norm(c,1) ) 
                pass(j1,j2,j3) = 0; 
            end
            j1 = j1 + 1; 
        end
        j2 = j2 + 1; 
    end
    j3 = j3 + 1; 
end

if ( all(all(all( pass ) ) ) ) 
    pass = all(all(all(pass))); 
end

end