function [v,lambda,k] = rayleighQuotientIteration( A, v, make_plots )
% rayleighQuotientIteration( A )
% Given a mxm real symmetric matrix A, returns the largest eigenvalue
% of A along with the corresponding eigenvector
%

% get the dimensions of A
[m,~] = size( A );

if nargin < 3
   make_plots = false; 
end

if nargin < 2
    
    % create a random initial vector
	v = rand( m, 1 );

    % normalize the initial vector
    v = v / norm( v );

end

% calculate initial lambda guess using Rayleigh Quotient
lambda = v'*A*v;

% the maximum number of iterations to perform
max_iterations = 100;

% algorithm will stop when the difference between
% successive eigenvalues falls below the value
stop_threshold = 1e-8;

% store all the previous eigenvalues, for plotting
% and for stopping criteria calculation
lambdas = zeros( max_iterations, 1 );

for k = 1:max_iterations

    w = ( A - lambda*eye(m) ) \ v;
    v = w / norm( w );
    lambda = v'*A*v;
    lambdas(k) = lambda;
    
    % check the stopping condition
    if ( k > 1 )
        d_lambda = abs( lambdas(k)-lambdas(k-1) );
        if ( d_lambda < stop_threshold )
            break;
        end
    end
    
end

lambda = lambdas(k);

if ( make_plots )
    plotResults( A, k, lambdas );
end

end
