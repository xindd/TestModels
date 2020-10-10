function [a,b,c,Error,Error2] = WWWCost(N,tL,PR,TR,PL,TL,TimeGrid)

if nargin < 7
    TimeGrid = [0,2];
    Length = 2;
end
% tL = 1/6
% N = 30
% TimeGrid = [2, 4]
% ±¨´í
% PL = [9.591286 9.370115];
% TL = [15.58621 16.00345];
% PR = [9.196904 9.443084];
% TR = [29.63946 35.54582];

% ²¨¶¯´ó
% PL = [0.4488358 0.5062905];
% TL = [1.531284 1.755409];
% PR = [0.7655, 0.7952];
% TR = [2.749777 2.334962];
if nargin < 3
    PL = [0.1869, 0.3816];
    TL = [0.4898, 0.6463];
    PR = [0.4387, 0.4456];
    TR = [0.7655, 0.7952];
end

if nargin < 2 
    tL = 1;
end

if nargin < 1
    N = 3;
end

Stop1 = 1.0;
Stop2 = 1.0;

PError = 0;
PError2 = 0;

alpha = 10000*ones(2,1);

a = zeros(1,N+1);
b = zeros(1,N+1);
c = ones(1,N+1);
h = (TimeGrid(end) - TimeGrid(1))/N;
Time = h : h : TimeGrid(end)-h;
a(1) = TL(1)/tL; a(end) = TL(end)/tL;

a = a(1) :(a(end) - a(1))/N : a(end);

Coef = [tL*PL(1)/TL(1),tL*PL(end)/TL(end)];

f = @(x) (1 - exp(-tL*x))./x - Coef;
g = @(x) (tL.*x.*exp(-tL*x) - (1 - exp(- tL*x)))./x.^2;
Temp =  [c(1),c(end)];
for i = 1 : 10
    Temp = Temp - (f(Temp))./g(Temp);
end

c = Temp(1) :(Temp(end) - Temp(1))/N : Temp(end);

IA = diag(2*ones(1,N-1)) + diag(-1 * ones(1,N-2),1)+ diag(-1 * ones(1,N-2),-1);
IB = diag(2*ones(1,N+1)) + diag(-1 * ones(1,N),1)+ diag(-1 * ones(1,N),-1);
IB(1,1) = 1; IB(end,end) = 1;
IC = IA;

MatrixJ = blkdiag(IA,IB,IC);
VectorJ = zeros(3*N-1,1);
VectorJ(1) = -a(1); VectorJ(N-2) = -a(N+1); VectorJ(2*N+1) = -c(1); VectorJ(end) = -c(end);
X = [a(2:end-1),b,c(2:end-1)];

for i = 1 : 50
    i=1;
    Grad = MatrixJ * X' + VectorJ;
    [S,GRAD] = PRTRX(a,b,c,N,TimeGrid(1:2),[PR(1);TR(1)]);
    Grad = Grad + alpha(1) * (S(1) - PR(end)) * GRAD(:,1) + alpha(2) * (S(2) - TR(end)) * GRAD(:,2);
    Index = find(Grad>0);
    rho = 0.5 * min(X(Index)./Grad(Index)');

    X = X - rho * Grad';
    StopRule = min([Stop1,Stop2,norm(Grad)]);
    if  StopRule < 1e-4
        break;
    end
    a(2:end-1) = X(1:N-1);
    b = X(N:2*N);
    c(2:end-1) = X(2*N+1:3*N-1);

    Error = X * MatrixJ * X'/2/h + X * VectorJ/h;
    Error2 = (S(1) - PR(end))^2 + (S(2) - TR(end))^2;
    Stop1 = abs(Error-PError);
    Stop2 = abs(Error2-PError2);
    PError = Error;
    PError2= Error2;

end
Error = X * MatrixJ * X'/2/h + X * VectorJ/h;
Error2 = (S(1) - PR(end))^2 + (S(2) - TR(end))^2;

plot(Time,a(2:end-1),'*',Time,b(2:end-1),'.',Time,c(2:end-1));
legend('a','b','c')
hold off

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