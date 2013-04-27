function [Q] = formQ( W )
% formQ( W ) Given a the W matrix produced by house( A ), build the
%            mxm orthogonal matrix Q which, together with the upper
%            triangular R matrix produced by house( A ), forms the
%            QR factorization of A
%


% get the dimensions of W
[m,n] = size(W);

Q = zeros(m,m);

% the basis vectors we will use to compute Q one vector at a time
% (essentially computing Q * I )
X = eye(m);

% outer loop computes Q*e1, Q*e2, ...
for i = 1:m
    
    x = X(:,i);
    
    % inner loop applies calculated Householder reflection vectors
    for k = n:-1:1
        
        % pull Householder reflection vectors out of W
        v = W(k:m,k);
        
        % apply v to x
        x(k:m) = x(k:m) - 2*(v*v')*x(k:m);
        
    end
    
    % store the calculated column of Q
    Q(:,i) = x;
    
end

end
