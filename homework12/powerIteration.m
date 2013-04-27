function [v,lambda] = powerIteration( A )
% powerIteration( A )
% Given a mxm real symmetric matrix A, returns the largest eigenvalue
% of A along with the corresponding eigenvector
%

% get the dimensions of A
[m,~] = size( A );

% create a random initial vector
v = rand( m, 1 );

% normalize the initial vector
v = v / norm( v );

% the maximum number of iterations to perform
max_iterations = 100;

% algorithm will stop when the difference between
% successive eigenvalues falls below the value
stop_threshold = 1e-10;

% store all the previous eigenvalues, for plotting
% and for stopping criteria calculation
lambdas = zeros( max_iterations, 1 );

for k = 1:max_iterations

    w = A * v;
    v = w / norm( w );
    lambdas(k) = v'*A*v;
    
    % check the stopping condition
    if ( k > 1 )
        d_lambda = abs( lambdas(k)-lambdas(k-1) );
        if ( d_lambda < stop_threshold )
            break;
        end
    end
    
end

lambda = lambdas(k);

% calculate the true result for plotting purposes
[~,D] = eig( A );
true_lambda = D(m,m);

% plot the error between the approximation of the eigenvalue
% and the "true" eigenvalue returned by matlab
figure;
semilogy( 1:k, abs(lambdas(1:k)-true_lambda),...
                '--rs','LineWidth',2,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','g',...
                'MarkerSize',10 );
title( 'Power Iteration: Error in Calculated Eigenvalue by Iteration' );
ylabel( 'Error' );
xlabel( 'Iteration' );
            
end