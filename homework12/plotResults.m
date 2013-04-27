function plotResults( A, k, lambdas )
% Helper function for plotting results generated by powerIteration.m
% and rayleighQuotientIteration.m

% get the dimensions of A
[m,~] = size( A );

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
title( 'Rayleigh Quotient Iteration: Error in Calculated Eigenvalue by Iteration' );
ylabel( 'Error' );
xlabel( 'Iteration' );

% plot the lambda ratios, which approximate the square of the
% true ratio of the second and first largest eigenvalues

lambda_ratio = ( lambdas(3:k) - lambdas(2:k-1) ) ./ ...
               ( lambdas(2:k-1) - lambdas(1:k-2) )

true_ratio = D(m-1,m-1)/D(m,m)

figure;
hold on;
plot( 3:k, lambda_ratio,...
                '--rs','LineWidth',2,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','g',...
                'MarkerSize',10 );
title( 'Rayleigh Quotient Iteration: Lambda Ratio by Iteration' );
ylabel( 'Ratio' );
xlabel( 'Iteration' );
plot([3,k],[true_ratio^2,true_ratio^2],...
                'LineWidth',1,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','g');
hold off;
end