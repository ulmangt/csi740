%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Experiment 1                  %%%
%%% Replicate the original result %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%randn('seed',232)

disp 'Experiment 1'

R = triu( randn(50) );
[Q,X] = qr( randn(50) );
A = Q*R;
[Q2,R2] = qr( A );

% poor accuracy
norm(Q2-Q)
norm(R2-R)/norm(R)

% accuracy near machine epsilon
norm(A-Q2*R2)/norm(A)

% construct another approximation of Q and R which
% individually are better than the Q2 and R2 produced
% by householder triangularization
Q3 = Q + (1e-4)*randn(50);
R3 = R + (1e-4)*randn(50);

% but the error is large now
norm(A-Q3*R3)/norm(A)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Experiment 2                  %%%
%%% Ensure Q3 is unitary and R3   %%%
%%% is upper triangular           %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp 'Experiment 2'

% try ensuring Q3 and R3 are orthogonal / upper triangular
% another qr factorization on Q3 produces a similar, but unitary matrix
[Q4 X] = qr( Q3 );
% simply zero out the lower diagonal elements of R3
R4 = triu(ones(50)).*R3; 

% the error should still be large
norm(A-Q4*R4)/norm(A)

% the errors in Q2 and R2 are forward errors resulting from the high
% condition number of R2 (random triangular matrices are ill conditioned
% with high probability
cond(R2)

% the error in Q2*R2 is a backwards error, suggesting the householder
% triangularization algorithm is backwards stable

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Experiment 3                  %%%
%%% Use a better conditioned R    %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp 'Experiment 3'

R = triu( ones(50) );
[Q,X] = qr( randn(50) );
A = Q*R;
[Q2,R2] = qr( A );

cond( R )

% accuracy near machine epsilon
norm(Q2-Q)
norm(R2-R)/norm(R)

% accuracy near machine epsilon
norm(A-Q2*R2)/norm(A)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Experiment 4                                              %%%
%%% Perform the experiment using single precision R and Q     %%%
%%% but calculate A, R2, and Q2 "exactly" in double precision %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp 'Experiment 4'

R = single( triu( randn(50) ) );
[Q,X] = qr( rand(50) );
Q = single( Q );
A = double( Q ) * double( R );
[Q2,R2] = qr( A );

% poor accuracy
norm(Q2-Q)
norm(R2-R)/norm(R)

% accuracy near machine epsilon
norm(A-Q2*R2)/norm(A)

% here we see that even knowing the exact QR factorization of A does
% not reduce the forward error

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Experiment 5                  %%%
%%%                               %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp 'Experiment 5'

R = [ 1 2 ; 0 3 ];
theta = -0.1;
Q = [ cos(theta) -sin(theta) ; sin(theta) cos(theta) ];
A = Q * R;
[Q2,R2] = qr( A );

Q2(:,1)*=-1;
R2(1,:)*=-1;

norm(Q2-Q)
norm(R2-R)/norm(R)
norm(A-Q2*R2)/norm(A)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Experiment 6                  %%%
%%%                               %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp 'Experiment 6'

R = [ 1 1 ; 0 eps ];
theta = -0.1;
Q = [ cos(theta) -sin(theta) ; sin(theta) cos(theta) ];
A = Q * R;
[Q2,R2] = qr( A );

norm(Q2-Q)
norm(R2-R)/norm(R)
norm(A-Q2*R2)/norm(A)
