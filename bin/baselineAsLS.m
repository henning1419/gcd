function z = baselineAsLS(y, lambda, p)
% Estimate baseline with asymmetric least squares
y=y(:);
m = length(y);
D = diff(speye(m), 2);
w = ones(m, 1);
for it = 1:10
    W = spdiags(w, 0, m, m);
    C = chol(W + lambda * (D' * D));
    z = C \ (C' \ (w .* y));
    w = p * (y > z) + (1 - p) * (y < z);
end
end