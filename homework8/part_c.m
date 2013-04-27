
format long

% the norm to use
n = 1;

% the b matrix from the problem statement
b = [ 1 ; 0 ];

% the number of epsilon values to investigate
steps = 16;

% allocate space for the calculated values
errors = zeros( steps, 1 );
relative_errors = zeros( steps, 1 );
norm_calc_x = zeros( steps, 1 );
norm_exact_x = zeros( steps, 1 );
condition_numbers_exact = zeros( steps, 1 );
condition_numbers_calc = zeros( steps, 1 );
epsilons = zeros( steps, 1 );
exponents = 0:steps;

% start at machine epsilon and increase by factors of two until reaching 1
% as expected, because of the 52 bits in the mantissa of an IEEE double,
% this takes 52 steps
for exponent=exponents
    
    % calculate epsilon
    %epsilon=eps*10.^exponent;
    epsilon=2*10^-exponent;
    
    % form the A matrix
    A = [ 1 1+epsilon ; 1-2*epsilon 1 ];
    
    % solve the equation
    x = A \ b;
    
    % compute the exact solution
    denominator = 2 * epsilon^2 + epsilon ;
    x_exact = [ 1  ;  2 * epsilon - 1 ] ./ denominator;
    
    % calculate the solution relative error
    
    error = norm(x-x_exact,n);
    %relative_error = norm(x,n) / norm(x_exact,n) - 1;
    relative_error = error / norm(x_exact,n);
    
    norm_exact_x(exponent+1) = norm(x_exact,n);
    norm_calc_x(exponent+1) = norm(x,n);
    epsilons(exponent+1) = epsilon;
    errors(exponent+1) = error;
    relative_errors(exponent+1) = relative_error;
    condition_numbers_exact(exponent+1) = ( 2 + epsilon )^2 / denominator;
    condition_numbers_calc(exponent+1) = cond( A, n );
end

loglog( epsilons, relative_errors );
title( 'Relative Error of x as a function of Epsilon' );
xlabel( 'epsilon' );
ylabel( 'Relative Error' );

% figure;
% semilogx( epsilons, errors );
% title( 'Absolute Error of x as a function of Epsilon' );
% xlabel( 'epsilon' );
% ylabel( 'Absolute Error' );
% 
% figure;
% loglog( epsilons, norm_calc_x ./ norm_exact_x );
% title( 'norm( calculated x ) / norm( exact x )' );
% xlabel( 'epsilon' );
% ylabel( 'Norm Ratio' );
% 
% figure;
% loglog( epsilons, condition_numbers_exact );
% title( 'Condition Number of A as a function of Epsilon' );
% xlabel( 'epsilon' );
% ylabel( 'Condition Number' );
% 
% figure;
% loglog( epsilons, abs( condition_numbers_exact - condition_numbers_calc ) );
% title( 'Difference between calcualted and exact Condition Number of A' );
% xlabel( 'epsilon' );
% ylabel( 'Condition Number Difference' );
% 
% figure;
% loglog( epsilons, norm_exact_x );
% title( '2 Norm of x as a function of Epsilon' );
% xlabel( 'epsilon' );
% ylabel( 'norm(x,2)' );