

function [b,LL] = fast_regress(y,X,bonly)

% function [b,LL] = fast_regress(y,X,bonly)
%
% standard regression function reduced to essentials to save processing
% time.
% returns beta and LL if needed. 

LL = NaN;

[n,ncolX] = size(X);

% Use the rank-revealing QR to remove dependent columns of X.
[Q,R,perm] = qr(X,0);
p = sum(abs(diag(R)) > max(n,ncolX)*eps(R(1)));
if p < ncolX
    R = R(1:p,1:p);
    Q = Q(:,1:p);
    perm = perm(1:p);
end

% Compute the LS coefficients, filling in zeros in elements corresponding
% to rows of X that were thrown out.
b = zeros(ncolX,1);
b(perm) = R \ (Q'*y);

if bonly % if we don't need the LL
  return;
end

yhat = X*b;                     % Predicted responses at each data point.
r = y-yhat;                     % Residuals.
sigma = std(r);

% log likelihood
prod = (yhat-y)'*(yhat-y);
LL = -n*log(2*pi)/2 - n*log(sigma^2)/2 - (1/(2*sigma.^2))*prod;


return;
