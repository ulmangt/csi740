% problem 10.3

Z = [ 1 2 3 ;
      4 5 6 ;
      7 8 7 ;
      4 2 3 ;
      4 2 2 ];

% get the dimensions of W
[m,n] = size(Z);
  
% modified Gram-Schmidt computation
[mgs_Q, mgs_R] = mgs( Z );

% MATLAB computation
[qr_Q, qr_R] = qr( Z, 0 );

% Householder triangularization computation
[house_W, house_R] = house( Z );
house_Q = formQ( house_W );
% get reduced Q
house_Q = house_Q(:,1:n);

% create a function to compute the mean of the absolute values of the
% differences between the components of two matrices
diff = @(A) (mean(abs(A(:))));

% qr / mgs difference
fprintf( 'MSG / QR  Q Error: %0.5e\n', diff( mgs_Q - qr_Q ) );
fprintf( 'MSG / QR  R Error: %0.5e\n', diff( mgs_R - qr_R ) );

% qr / house difference
fprintf( 'QR / House  Q Error: %0.5e\n', diff( house_Q - qr_Q ) );
fprintf( 'QR / House  R Error: %0.5e\n', diff( house_R - qr_R ) );

% house / mgs difference
fprintf( 'House / MGS  Q Error: %0.5e\n', diff( house_Q - mgs_Q ) );
fprintf( 'House / MGS  R Error: %0.5e\n', diff( house_R - mgs_R ) );


fprintf( 'House Reconstructed A Error: %0.5e\n', diff( house_Q*house_R - Z ) );
fprintf( 'MGS Reconstructed A Error: %0.5e\n', diff( mgs_Q*mgs_R - Z ) );
fprintf( 'QR Reconstructed A Error: %0.5e\n', diff( qr_Q*qr_R - Z ) );