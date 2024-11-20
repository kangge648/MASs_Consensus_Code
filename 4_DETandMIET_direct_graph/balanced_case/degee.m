clear
clc

% compute the guaranteed MIETs
A = [0.1 1
    0  0];
B_u = [-1;1];
ADJ  = [0 0 0 0 1
    1 0 0 0 0
    0 1 0 0 0
    1 0 1 0 0
    0 0 0 1 0];
L = diag(sum(ADJ,2))-ADJ;
H = diag([2 1 1 1 2])/7;
h = [2 1 1 1 2]'/7;

phi = 0.618;
psi = norm(H^0.5*L*H^-0.5);
xishu = psi^2*(1+psi/phi)^2;

N = length(L);
nx = size(A,1);
nu = size(B_u,2);

Q = eye(nx);
P = are(A, B_u*B_u', Q); % are(A,B,C): A'*X + X*A - X*B*X + C = 0
K = -B_u'*P;

lambda = eig(K'*K/P);
lambda_k = max(lambda);
lambda = eig(Q/P);
lambda_q = min(lambda(lambda>10^-10));
beta = 50;
delta = 0.01*[1 1 1 1 1];
for k = 1:20
    omega = 0.1*k;
    c = 1/phi/(1-omega);
    for i = 1:(N-1)
        Deta(k,i) = 4*lambda_k^2*xishu*c^2/omega-(delta(i)+lambda_k-lambda_q)^2;
        L_T(k,i) = 2/sqrt(Deta(k,i))*(atan((2*beta+delta(i)+lambda_k-lambda_q)/sqrt(Deta(k,i)))-atan((delta(i)+lambda_k-lambda_q)/sqrt(Deta(k,i))));
    end

    for i = N:N
        Deta(k,i) = 4*lambda_k^2*xishu*c^2/omega-(delta(i)+lambda_k-lambda_q)^2
        L_T(k,i) = 2/sqrt(Deta(k,i))*(atan((2*beta+delta(i)+lambda_k-lambda_q)/sqrt(Deta(k,i)))-atan((delta(i)+lambda_k-lambda_q)/sqrt(Deta(k,i))))
    end
end