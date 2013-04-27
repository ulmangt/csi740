function [W,R] = house( A )
% mgs( A ) Given a mxn complex matrix A, decomposes the matrix into
%          an nxn upper triangular matrix R and a lower
%          triangular matrix W consisting of n successive Householder
%          reflection vectors
%

% get the dimensions of A
[m,n] = size(A);

% initialize W and R matrices
W = zeros( m, n );

for k = 1:n

    % x is the vector which, when reflected over the Householder
    % reflection vector, results in a vector of equal length with
    % all zeros expect in the first index
    x = A(k:m,k);
    
    % build the Householder reflection vector v
    v = x;
    v(1) = v(1) + sign(x(1))*norm(x);
    
    if ( norm(v) ~= 0 )
        v = v / norm(v);
    end
    
    % and store it in the lower triangular matrix W
    W(k:m,k) = v;
    
    % update the kth column of A using v
    A(k:m,k:n) = A(k:m,k:n) - 2*(v*v')*A(k:m,k:n);
    
end

% return the square and upper triangular portion of A
R = A(1:n,1:n);

end