function [A,errors_fro,errors_max,k] = qrIteration( A )

% get the dimensions of A
[m,~] = size( A );

% build a mask selecting the lower triangular part of A
% this is used in the error calculations
mask = ~triu( ones(m) );

% the maximum number of iterations to perform
max_iterations = 100;

% algorithm will stop when the difference between
% successive eigenvalues falls below the value
stop_threshold = 1e-12;

% store all the previous eigenvalues, for plotting
% and for stopping criteria calculation
errors_fro = zeros( max_iterations, 1 );
errors_max = zeros( max_iterations, 1 );

for k = 1:max_iterations

    % perform one QR iteration step
    [Q,R] = qr( A );
    A = R*Q;
    
    % multiply A element-wise by the pre-calculated mask
    % then take the frobenius norm of the remaining matrix
    errors_fro(k) = norm( mask .* A, 'fro' );
    
    % also estimate the error by simply taking the largest value
    errors_max(k) = max(max( mask .* A ));
     
    % check the stopping condition
    if ( errors_max(k) < stop_threshold )
       break 
    end
    
end

end