
%%%%%%%%% Problem 5

%%%%%%%%% Part a

% input data points
t = [ 4 6 7 8 9 13 14 15 16 17 18 19 20 21 26 ]';
y = [ 1319 1197 1075 1101 1010 675 715 643 581 534 513 480 446 443 328 ]';

% plot data points
scatter( t, y, '+' ); hold on;

% build Vandermond matrix
A = [ t.^0 t.^1 ];

% compute QR factorization of Vandermon matrix A
[Q R] = qr( A );

% solve the upper triangular system Rc = Q'y for c
c = R \ (Q'*y);

% c holds the coefficients of the linear fit y = c(2)*x + c(1)
% print the linear fit equation
fprintf( '\ny = %g x + %g\n', c(2), c(1) );

% find the x-intercept
t0 = -c(1)/c(2);
fprintf( '\nx=%g y=0\n', t0 );

% plot the linear fit extended to the x and y intercepts
x = [ 0 ; t ; t0 ];
plot( x, c(2)*x+c(1) );

%%%%%%%%% Part b

% plot data points
figure; scatter( t, y, '+' ); hold on;

% compute SVD factorization of A
[U,S,V] = svd( A );

% solve the diagonal system S*w = U'*y
w = S \ (U'*y);

% compute the coefficient vector c
c = V*w;

% c holds the coefficients of the linear fit y = c(2)*x + c(1)
% print the linear fit equation
fprintf( '\ny = %g x + %g\n', c(2), c(1) );

% find the x-intercept
t0 = -c(1)/c(2);
fprintf( '\nx=%g y=0\n', t0 );

% plot the linear fit extended to the x and y intercepts
x = [ 0 ; t ; t0 ];
plot( x, c(2)*x+c(1) );