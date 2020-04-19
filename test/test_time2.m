% execution time of different approaches to obtain eigenvalues(Iteration time: 500)
function test_time2()

addpath(genpath('../src/'))
s = zeros(1,10);
t1 = zeros(1,10);
t2 = zeros(1,10);
t3 = zeros(1,10);
n_iter = 500;
for i=1:10
  A = rand(i*5);
  tic;
  pure_qr(A, n_iter, "hh");
  t1(i) = toc;
  % Placeholder for practical qr
  tic;
  eig(A);
  t2(i) = toc;
  
  tic;
  eig(A);
  t3(i) = toc;

  s(i)= i*5;
end

figure
ax1 = subplot(1,3,1);
p1 = plot(s,t1,'LineWidth',1);
ax1.XGrid = 'on';
ax1.YGrid = 'on';
p1.Marker = 'o';
title('Pure QR');
ylabel('Time (sec)');
xlabel('Matrix size (nxn)');

ax2 = subplot(1,3,2);
p2 = plot(s,t2,'r','LineWidth',1);
ax2.XGrid = 'on';
ax2.YGrid = 'on';
p2.Marker = 'o';
title('Practial QR');
ylabel('Time (sec)');
xlabel('Matrix size (nxn)');

ax3 = subplot(1,3,3);
p3 = plot(s,t3,'g','LineWidth',1);
ax3.XGrid = 'on';
ax3.YGrid = 'on';
p3.Marker = 'o';
title('Built-in Eig function');
ylabel('Time (sec)');
xlabel('Matrix size (nxn)');
print('Runtimes Difference across algorithms for obtaining eigenvalues', '-depsc');
 

end