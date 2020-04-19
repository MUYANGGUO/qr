% Comparison of orthogonality of Q matrix for different QR decomp approaches(Fixed: QR algorithm: Pure QR; Fixed matrix size Iteration time: 500)
function test_orthogonality1()
n_iter = 500;
s = zeros(1,n_iter);
t1 = zeros(1,n_iter);
t2 = zeros(1,n_iter);
t3 = zeros(1,n_iter);
A0 = rand(50,50);
% [~,n] = size(A0);
I = eye(size(A0));
format long
% 
%% HH
A = A0;
for i=1:n_iter
    [Q,R] = houseqr(A);
    ortherr_hh = norm(Q'*Q-I,inf);
    t1(i) = ortherr_hh;
	A = R * Q;
    s(i) = i;
end
%% Givens
A = A0;
for i=1:n_iter
    [Q,R] = givensqr(A);
    ortherr_givens = norm(Q'*Q-I,inf);
    t2(i) = ortherr_givens;
	A = R * Q;
end

%% mgs
A = A0;
for i=1:n_iter
    [Q,R] = mgsqr(A);
    ortherr_mgs = norm(Q'*Q-I,inf);
    t3(i) = ortherr_mgs;
	A = R * Q;
end
%%
disp(ortherr_hh);
disp(ortherr_givens);
disp(ortherr_mgs);


figure
ax1 = subplot(1,3,1);
p1 = plot(s,t1,'LineWidth',1);
ax1.XGrid = 'on';
ax1.YGrid = 'on';
p1.Marker = 'o';
title('Householder QR');
ylabel('Othogonality Error');
xlabel('Number of Iterations');

ax2 = subplot(1,3,2);
p2 = plot(s,t2,'r','LineWidth',1);
ax2.XGrid = 'on';
ax2.YGrid = 'on';
p2.Marker = 'o';
title('Givens QR');
ylabel('Othogonality Error');
xlabel('Number of Iterations');

ax3 = subplot(1,3,3);
p3 = plot(s,t3,'g','LineWidth',1);
ax3.XGrid = 'on';
ax3.YGrid = 'on';
p3.Marker = 'o';
title('Modified Gram-Schmit');
ylabel('Othogonality Error');
xlabel('Number of Iterations');
print('Othogonality Error of Different QR Decomposition', '-depsc');

end

