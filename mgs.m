function [Q,R] = mgs( A )
% mgs( A ) Given a mxn complex matrix A, decomposes the matrix into
%          a mxn matrix with orthonormal columns Q and an nxn upper
%          upper triangular matrix R using the modified Gram-Schmidt
%          algorithm
%

% get the dimensions of A
[m,n] = size(A);

% initialize V, Q, and R matrices
V = A;
Q = zeros( m, n );
R = zeros( n, n );

% iterate over the columns of A, building columns of Q
for i = 1:n
    
    % build ith diagonal of R and ith column of Q
    R(i,i) = norm(V(:,i));
    
    % if R(i,i) == 0 we simply pick any normalized vector orthogonal to
    % Q(:,1) through Q(:,i-1) and continue
    if ( R(i,i) == 0 )
        
        Q(:,i) = choose_ortho_vector( Q, i );
        
    else
    
        Q(:,i) = V(:,i) / R(i,i);
    
    end
    
    % remove components of remaining i+1 through n columns of V
    % in the direction of the ith column of Q
    for k = i+1:n
       
        R(i,k) = Q(:,i)'*V(:,k);
        V(:,k) = V(:,k) - R(i,k)*Q(:,i);
        
    end
    
end

end

function basis = choose_ortho_vector( Q, i )
% choose_ortho_vector( Q, i ) a helper function which deals with the case
%                             R(i,i) == 0 by choosing an arbitrary
%                             normalized vector orthogonal to the existing
%                             Q vectors
%

    [m,~] = size(Q);

    % create n basis vectors
    all_basis = eye(m);

    % special case i == 1
    if ( i == 1 )

        basis = all_basis(:,1);

    else

        % check each basis vector against Q(:,1) through Q(:,i-1)
        for j = 1:m

            if ( zeros(i-1,1) == Q(:,1:i-1)'*all_basis(:,j) )
                break
            end

        end

        basis = all_basis(:,j);

    end
end