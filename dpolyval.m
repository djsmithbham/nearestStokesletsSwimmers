function dy=dpolyval(p,x)

% evaluates first derivative of polynomial with coefficients in descending
% order in p, at x
p=p(:)';
q=p(1:end-1).*[length(p)-1:-1:1];
dy=polyval(q,x);
