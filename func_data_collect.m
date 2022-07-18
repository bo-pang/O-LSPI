function [zz,zx,zt_r,A,B,C,Q,R,K...
    ] = func_data_collect(T,noise_magnitude)

rng shuffle;

A = [0.95 0.01 0;
    0.01 0.95 0.01;
    0 0.01 0.95];
B = [1 0.1;0 0.1;0 0.1;];
[n,m] = size(B);
C = eye(n);

Q = eye(n);
R = eye(m);

n1 = m+n+1;
n2 = n1*(n1+1)/2;
n3 = n*(n+1)/2;

% [K,~] = place(A,B,[0.8,0.9]);
K = zeros(m,n);
x = zeros(n,1);
u = zeros(m,1);
zz = zeros(n2,n2);
zx = zeros(n2,n3);
zt_r = zeros(n2,1);
for i=1:T
    w = randn(m+n,1);
    u = -K*x + noise_magnitude*w(n+1:end);
    zt = kronv([x;u;1]);
    zz = zz + zt*zt';
    x_next = A*x + B*u + C*w(1:n); 
    zx = zx + zt*kronv(x_next)';
    zt_r = zt_r + zt*(x'*Q*x+u'*R*u);
    x = x_next;
end
end