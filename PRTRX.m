function [S,GRAD] = PRTRX(a,b,c,N,TimeGrid,Initial)

eta = zeros(3*N-1,N);
for i = 2 : N    
    eta(i-1,i-1) = 1;
end

xi = zeros(3*N-1,N);
for i = 2 : N        
    xi(2*N+i-1,i-1) = 1;
end

theta = zeros(3*N-1,N);
for i = 1 : N        
    theta(N-1+i,i) = 1;
end

f = @(x,a,b,c) [a - c * x(1); a - b * (x(2)-x(1))];
g = @(x,eta,xi,theta,b,c,S) [ eta - xi * S(1) - c * x(1:3*N-1) ; 
                              eta - (S(2)-S(1))*theta - b*(x(3*N:end)-x(1:3*N-1))];

h = (TimeGrid(end) - TimeGrid(1))/N;
S = Initial;
K = zeros(2,4);
GRAD = zeros(6*N-2,1);
L = zeros(6*N-2,4);

for i = 1 : N - 1
    
    S1 = S;
    K(:,1) = f(S1,a(i),b(i),c(i));
    S2 = S + h/2 * K(:,1);
    K(:,2) = f(S2,a(i),b(i),c(i));
    S3 = S + h/2 * K(:,2);
    K(:,3) = f(S3,a(i),b(i),c(i));
    S4 = S + h * K(:,3);
    K(:,4) = f(S4,a(i),b(i),c(i));
    S = S + h * ( K(:,1) + 2 * K(:,2) + 2 * K(:,3) + K(:,4))/6;
    
    L(:,1) = g(GRAD,eta(:,i),xi(:,i),theta(:,i),b(i),c(i),S1);
    L(:,2) = g(GRAD+h/2 * L(:,1),eta(:,i),xi(:,i),theta(:,i),b(i),c(i),S2);
    L(:,3) = g(GRAD+h/2 * L(:,2),eta(:,i),xi(:,i),theta(:,i),b(i),c(i),S3);
    L(:,4) = g(GRAD+h * L(:,3),eta(:,i),xi(:,i),theta(:,i),b(i),c(i),S4);
    GRAD = GRAD + h * ( L(:,1) + 2 * L(:,2) + 2 * L(:,3) + L(:,4))/6;
    
end

GRAD = reshape(GRAD,3*N-1,2);