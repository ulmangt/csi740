
% Matlab driver code to apply the qrIteration command to the provided data

% input the A matrix
A = [ 1 1  1   1   1   1 ;
      1 2  3   4   5   6 ;
      1 3  6  10  15  21 ;
      1 4 10  20  35  56 ;
      1 5 15  35  70 126 ;
      1 6 21  56 126 252  ];

[A_prime1, errors_fro, errors_max, k] = qrIteration( A );

figure;
semilogy( 1:k, errors_max(1:k),...
                '--rs','LineWidth',2,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','g',...
                'MarkerSize',6 );
title( 'QR Iteration' );
ylabel( 'Size of Largest Lower Triangular Element' );
xlabel( 'Iteration' );

% make the A matrix no longer symmetric
A(1,6) = 0;

[A_prime2, errors_fro, errors_max, k] = qrIteration( A );

figure;
semilogy( 1:k, errors_max(1:k),...
                '--rs','LineWidth',2,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','g',...
                'MarkerSize',6 );
title( 'QR Iteration on non-symmetric A' );
ylabel( 'Size of Largest Lower Triangular Element' );
xlabel( 'Iteration' );