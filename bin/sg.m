%% sg.m - Savitzky-Golay length-N order-d smoother design.
%
% [B, S] = sg(d, N);
%
% N = 2 M + 1 = filter length, d = polynomial order
% S = [ s 0 , s 1 , . . . , s d ] , F = S T S
% G = SF − 1 = derivative filters
% B = SF − 1 S T = smoothing filters
% indexing: B(M + 1 + m, M + 1 + k) = B mk , m, k = −M : M
% m-th SG filter = B( : , M + 1 + m) = b m , m = −M : M
% NRRs = diagonal entries of B .

function [B, G, S] = sg(d, N) %d übergeben wenn mit doppelter for-Schleife
M = (N-1)/2;
S = zeros(N,d+1);
for m=-M:M
    for i=0:d
        S(m+M+1, i+1) = m^i;
    end
end

% fliplr(vande(x)) laut Run&Time Analyse langsamer zur Erstellung von S
% v = -M:1:M;
% S2 = fliplr(vander(v));
% S2 = S2(:,1:d+1);
% 
% F2 = S2' * S2;
% G2 = S2 * F2^(-1);
% B2 = G2 * S2'; 

F = S' * S;
G = S * F^(-1);
B = G * S';
end

% polynomial coefficients c_i: c = G'*x 
% where G = S * F^(-1) = S * (S' * S)^(-1)