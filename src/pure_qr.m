function eig = pure_qr(A, n_iter, method)
%%
%Arguments:
%    n_iter: number of iteration for QR algorithm
%%
[m,n] = size(A);
eig = eye(n,1);

for i=1:n_iter
    if method == "hh"
        [Q,R] = houseqr(A);
    elseif method == "givens"
        [Q,R] = givensqr(A);
    else
        [Q,R] = mgsqr(A);
    end
	A = R * Q;
end

for i=1:n
	eig(i) = A(i,i);
end

%disp(A);

end