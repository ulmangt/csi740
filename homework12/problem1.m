
% input the A matrix
A = [ 1 1  1   1   1   1 ;
      1 2  3   4   5   6 ;
      1 3  6  10  15  21 ;
      1 4 10  20  35  56 ;
      1 5 15  35  70 126 ;
      1 6 21  56 126 252  ];

% perform power iteration to find the largest eigenvalue  
[v,lambda] = powerIteration( A )

% check the results using Matlab's built-in command
[V,D] = eig( A )