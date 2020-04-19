function compare()
% compare(X)
% Compare three QR decompositions,
%
X = rand(50);
I = eye(size(X));

% %% Classic Gram Schmidt
% [Q,R] = gs(X);
% qrerr_gs = norm(Q*R-X,inf)/norm(X,inf);
% ortherr_gs = norm(Q'*Q-I,inf);

%% Modified Gram Schmidt
[Q,R] = mgsqr(X);
qrerr_mgs = norm(Q*R-X,inf)/norm(X,inf);
ortherr_mgs = norm(Q'*Q-I,inf);

%% Householder QR Decomposition

[Q,R] = houseqr(X);
qrerr_hh = norm(Q*R-X,inf)/norm(X,inf);
ortherr_hh = norm(Q'*Q-I,inf);

% [U,R] = hqrd(X);
% QR = R;
% E = I;
% for k = size(X,2):-1:1
%     uk = U(:,k);
%     QR = QR - uk*(uk'*QR);
%     E = E - uk*(uk'*E) - (E*uk)*uk' + uk*(uk'*E*uk)*uk';
% end
% qrerr_h = norm(QR-X,inf)/norm(X,inf);
% ortherr_h = norm(E-I,inf);

%% Report results 
% fprintf('QR error\n')
% fprintf('Classic:     %10.3e\n',qrerr_gs)
% fprintf('Modified:    %10.3e\n',qrerr_mgs)
% fprintf('Householder: %10.3e\n',qrerr_h)
fprintf('\n')
fprintf('Orthogonality error\n')
% fprintf('Classic:     %10.3e\n',ortherr_gs)
% fprintf('Modified:    %10.3e\n',ortherr_mgs)
% fprintf('Householder: %10.3e\n',ortherr_hh)
disp(ortherr_mgs*10e14)
disp(ortherr_hh*10e14)