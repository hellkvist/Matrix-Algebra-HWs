clearvars
A = [1 2 3;
     2 4 5;
     3 5 6];
% find eigenvalues using eig
E = eig(A);
for i = 1:5
    A(:,:,i+1) = plane_rot_largest_off_diagonal(A(:,:,end));
end
disp(A(:,:,end))
disp(E)
function B = plane_rot_largest_off_diagonal(A)
    [n,~] = size(A);
    U = eye(n);
    
    [i,j] = idx_largest_off_diagonal(A);
    t = 0.5 * acot( (A(i,i) - A(j,j))/(2*A(i,j)) );
    
    ct = cos(t);
    st = sin(t);
    
    U(i,i) = ct;
    U(i,j) = -st;
    U(j,i) = st;
    U(j,j) = ct;
    
    B = U' * A * U;
end

function [i,j] = idx_largest_off_diagonal(A)
    % returns i,j where i < j, for the largest (magnitude)
    % off diagonal entry in A
    % A must be symmetric
    if all(all(A - A' < 1e-15))
        M = abs(A - diag(diag(A)));
        [a1, b1] = max(M);
        [~, j] = max(a1);
        i = b1(j);
        if i >= j
            [i,j] = deal(j,i);
        end
    else
        i = -1;
        j = -1;
    end
end