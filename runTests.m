function runTests( ) 
% RUNTESTS     Run all the tests
%
% Author:  
%    Alex Townsend, Jan 15 (originally written)

list = dir; 

for j = 1:length( list )
    if ( strfind(list(j).name,'test_') ) 
        eval( list(j).name(1:end-2) )
    end 
end

end