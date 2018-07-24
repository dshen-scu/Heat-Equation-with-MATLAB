% u(x, y, t)
% heat equation u_t = u_xx + u_yy + f
% initial data u(x, y, 0) = food court
% using finite-differences. Backward Time, Center Space.
% t_0 is the initial time, t_f is the final time
% k is the number of time steps
% m is the number of x space steps
% n is the number of y space steps
% !!!THIS CODE IS ONLY FOR THE CASE m=n!!!

function [u, x, y, t] = DTEST(t_0, t_f, k, m, n)
       p00 =        3581  ;
       p10 =      -22.55  ;
       p01 =      -12.89  ;
       p20 =       3.078  ;
       p11 =      0.2705  ;
       p02 =       2.988  ;
       p30 =     -0.2033  ;
       p21 =    -0.06079  ;
       p12 =     0.03021  ;
       p03 =     -0.3186  ;
       p40 =    0.006296  ;
       p31 =    0.003084  ;
       p22 =    -0.00139  ;
       p13 =  -3.353e-05  ;
       p04 =     0.01445  ;
       p50 =   -7.32e-05  ;
       p41 =  -5.234e-05  ;
       p32 =    4.21e-05  ;
       p23 =  -2.871e-05  ;
       p14 =   1.159e-05 ; 
       p05 =  -0.0002411;

       syms x y;
       f(x, y) = p00 + p10*x + p01*y + p20*x^2 + p11*x*y + p02*y^2 + p30*x^3 + p21*x^2*y + p12*x*y^2 + p03*y^3 + p40*x^4 + p31*x^3*y + p22*x^2*y^2 + p13*x*y^3 + p04*y^4 + p50*x^5 + p41*x^4*y + p32*x^3*y^2 + p23*x^2*y^3 + p14*x*y^4 + p05*y^5;
 
       
% Define Time Step
dt = (t_f-t_0)/k;
t = t_0:dt:t_f;

% Define x Space Step
dx = 28.8/(m-1);
x = 0:dx:28.8;
x = x';

% Define y Space Step
dy = 18.6/(n-1);
y = 0:dy:18.6;
y = y';

% Define ratio r
r = 0.25*dt/dx^2;

% Input Initial Data
u0 = zeros(m+2);


for i=1:m
    for j=1:n
        u0(i+1, j+1) = f(x(i), y(j));
    end
end

%Input Boundary Data
u0(1, 1) = f(x(1)-dx, y(1)-dy);
u0(1, n+2) = f(x(1)-dx, y(n)+dy);
u0(m+2, 1) = f(x(n)+dx, y(1)-dy);
u0(m+2, n+2) = f(x(n)+dx, y(n)+dy);
 
for i=1:m
    u0(i+1, 1) = f(x(i), y(1)-dy);
end
 
for j=1:n
    u0(1, j+1) = f(x(1)-dx, y(j));
end
 
for i=1:m
    u0(i+1, m+2) = f(x(i), y(n)+dy);
end
 
for j=1:n
    u0(m+2, j+1) = f(x(m)+dx, y(j));
end

% Define Matrix A
A = zeros((m)*(n));

a = 1+2*r; b= 1-2*r; c = -r/2;

P = zeros(m);
P(1, 1) = a; P(1, 2) = c; P(m, m-1) = c; P(m, m) = a;
for i=2:m-1
    P(i, i-1) = c; P(i, i) = a; P(i, i+1) = c;
end

Q = zeros(m);                 
for i=1:m
    Q(i, i) = c;
end

for i=1:n
    A((i-1)*m+1:i*m, (i-1)*m+1:i*m) = P;
end
for i=1:n-1
    A(i*m+1:(i+1)*m, (i-1)*m+1:i*m) = Q;
    A((i-1)*m+1:i*m, i*m+1:(i+1)*m) = Q;
end

% What to input on B (Ax=B)
B = zeros(m*n, 1);

% SOLVE!

% Answer Storage
u = zeros(m+2, m+2, k+1);
u0= 0.1*u0;
u(:, :, 1) = u0;

% External Force
F = zeros(m, m);

for i=1:m
    for j=1:n
        F(i, j) = 1*sin(0.7*x(i))*cos(0.7*y(j))*dt;
    end
end

T = zeros(m+2);


for l=1:k

%Boundary Force

U0 = 100*(1/1183.451)*u0-273;

for i=2:m+1
    T(i, 1) = abs(U0(i, 2)-U0(i, 1))*28.8*3.34*dt*3/(3600*m);
    T(i, m+2) = abs(U0(i, m+1)-U0(i, m+2))*28.8*3.34*dt*3/(3600*m);
end


for j=1:n
    for i=1:m
        B((j-1)*m+i, 1) = b*u0(i+1, j+1)-c*(u0(i+2, j+1)+u0(i, j+1)+u0(i+1, j+2)+u0(i+1, j));
    end
end

for j=1:n
    for i=1:m
        K((j-1)*m+i, 1) = F(i, j);
    end
end

% add K here
X = A^-1 * (B);
u1 = u0;

for j=1:n
    for i=1:m
    u1(i+1, j+1) = X((j-1)*m+i, 1);
    end
end

u(:, :, l+1) = u1;

u0 = u1+T;

end

for i=1:k+1
    surf(u(2:m+1, 2:m+1, i));
    axis([0 m+4 0 n+4 0 400]);
    pause(1);
end

end





