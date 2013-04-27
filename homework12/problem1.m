
% input the A matrix
A = [ 1 1  1   1   1   1 ;
      1 2  3   4   5   6 ;
      1 3  6  10  15  21 ;
      1 4 10  20  35  56 ;
      1 5 15  35  70 126 ;
      1 6 21  56 126 252  ];

% perform power iteration to find the largest eigenvalue  
[v,lambda,~] = powerIteration( A, rand(6,1), true )

% check the results using Matlab's built-in command
[V,D] = eig( A )

n_trials = 1000;
iteration_counts = zeros(n_trials,1);
for i = 1:n_trials
    
    [~,~,k] = powerIteration( A );
    iteration_counts(i) = k;
    
end

figure;
hist( iteration_counts, 3:1:8 );
title( 'Iterations Necessary For 1000 Random Starting Vectors' );
ylabel( 'Count' );
xlabel( 'Iterations' );