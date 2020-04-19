function [Q,R] = givensqr(A)
    [row,col] = size(A);
    Q = eye(row);
    R = A;
    for i = 1:col
        for j = i+1:row
            G = eye(row);
            c = R(i,i)/((R(i,i))^2 + (R(j,i)^2))^0.5;
            s = R(j,i)/((R(i,i))^2 + (R(j,i)^2))^0.5;
            G(i,i) = c;
            G(j,j) = c;
            G(j,i) = -s;
            G(i,j) = s;
            R = G*R;
            Q = G*Q;
        end
    end
    Q = Q';
end
