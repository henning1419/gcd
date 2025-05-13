%% sgfilt.m - filtering with length-N order-d SG smoother.
%
% y = sgfilt(d, N, x);
%
% x and y are L× 1 column vectors; and N = 2 M + 1. Must have L > N + 1.
% B( : , i) = b i−M− 1 = input-on transient filters, i = 1 : M + 1
% B( : , M + 1 ) = b 0 = steady-state filter
% B( : , M + 1 + i) = b i = input-off transient filters, i = 0 : M

function [y,y_abl1,y_abl2,y_abl3]= sgfilt(d, N, x)
M = (N-1)/2;
L = length(x);

% design filter
[B, G, ~]= sg(d, N);

% allocate memory
y = zeros(L,1);
y_abl1 = zeros(L,1);
y_abl2 = zeros(L,1);
y_abl3 = zeros(L,1);

%% input-on transients
for i = 1:M+1
    y(i,1) = B(:,i)' * x(1:N);
end 

%% steady state
b0 = B(:,M+1)';

X = zeros(N,L-2-2*M);
for n = 1:L-2*M-2
    X(:,n) = x(n+1:n+2*M+1);
end
y(M+2:L-M-1,1) = b0 * X;

% 1st derivative
if d >= 1
g1 = G(:,2)';
y_abl1(M+2:L-M-1,1) = g1 * X;
end

% 2nd derivative
if d >= 2
g2 = G(:,3)';
y_abl2(M+2:L-M-1,1) = 2*g2 * X;
end

% 3rd derivative
if d >= 3
g3 = G(:,4)';
y_abl3(M+2:L-M-1,1) = 6*g3 * X;
end

%% input-off transients
for i = 0:M
    y(L-M+i,1) = B(:,M+1+i)' * x(L-N+1:L);
end 

end
