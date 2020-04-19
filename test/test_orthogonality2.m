% Comparison of orthogonality of Q matrix for different QR decomp approaches(Fixed: QR algorithm: Pure QR; change matrix size? Fixed Iteration time: 500)
function test_orthogonality2()
n_iter = 500;
n_matrix = 10;
s = zeros(1,n_matrix);
t1 = zeros(1,n_matrix);
t2 = zeros(1,n_matrix);
t3 = zeros(1,n_matrix);

format long
% 

for k=1:n_matrix
  A0 = rand(k*5);
  s(k)= k*5;
  I = eye(size(A0));
  
  %% HH 
  A = A0;
  for i=1:n_iter
    [Q,R] = houseqr(A);
	A = R * Q;
  end
  ortherr_hh = norm(Q'*Q-I,inf);
  t1(k) = ortherr_hh;
  
  %% givens
  A = A0;
  for i=1:n_iter
    [Q,R] = givensqr(A);
	A = R * Q;
  end
  ortherr_givens = norm(Q'*Q-I,inf);
  t2(k) = ortherr_givens;
  
  %% mgs
  A = A0;
  for i=1:n_iter
    [Q,R] = mgsqr(A);
	A = R * Q;
  end
  ortherr_mgs = norm(Q'*Q-I,inf);
  t3(k) = ortherr_mgs;
end
  

figure
ax1 = subplot(1,3,1);
p1 = plot(s,t1,'LineWidth',1);
ax1.XGrid = 'on';
ax1.YGrid = 'on';
p1.Marker = 'o';
title('Householder QR');
ylabel('Othogonality Error');
xlabel('Matrix size (nxn)');

ax2 = subplot(1,3,2);
p2 = plot(s,t2,'r','LineWidth',1);
ax2.XGrid = 'on';
ax2.YGrid = 'on';
p2.Marker = 'o';
title('Givens QR');
ylabel('Othogonality Error');
xlabel('Matrix size (nxn)s');

ax3 = subplot(1,3,3);
p3 = plot(s,t3,'g','LineWidth',1);
ax3.XGrid = 'on';
ax3.YGrid = 'on';
p3.Marker = 'o';
title('Modified Gram-Schmit');
ylabel('Othogonality Error');
xlabel('Matrix size (nxn)');
print('Othogonality Error of Different Sizes of Matrix', '-depsc');

end

