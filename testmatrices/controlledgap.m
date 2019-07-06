function A = controlledgap(m,n,r,gap)
    %
    % Input
    % m,n - size of matrix
    % r   - location of gap
    % gap - size of gap
    
    % Output
    % A   - resulting matrix


    f = [gap./(1:r),1./(r+1:n)];

    A = sparse(m,n);
    for j = 1:n
        xj = sprand(m,1,0.025);
        yj = sprand(n,1,0.025);
        A = A + f(j)*xj*yj';
    end