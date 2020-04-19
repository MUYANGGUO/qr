% Return the Hessenberg form of the matrix
function hess = hessenberg(A)

[m,n] = size(A);
hess = A;

for i=1:n-2
	column = norm(hess(i+1:end,i),2) * eye(length(hess(i+1:end,i)),1) - hess(i+1:end,i);
	HH = eye(length(column)) - 2/(column' * column) * (column * column');
	B = eye(n);
	%disp(size(HH));
	%disp(size(B));
	B(i+1:end, i+1:end) = HH;
	hess = B * hess * B;
    disp(i)
end
end