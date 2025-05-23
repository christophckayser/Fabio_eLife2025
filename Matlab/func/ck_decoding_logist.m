function [v,lambda,LL] = ck_decoding_logist(x,y,v,show,regularize,lambda,lambdasearch,eigvalratio)

% logist() - Iterative recursive least squares algorithm for linear 
%	     logistic model - RETURNS Coefficients up to scale
%
% Usage:
%   >> [v] = logist(x,y,vinit,show,regularize,lambda,lambdasearch,eigvalratio);
% 
% Inputs:
%   x - N input samples [N,D]
%   y - N binary labels [N,1] {0,1}
%
% Optional parameters:
%   vinit       - initialization for faster convergence
%   show        - if>0 will show first two dimensions
%   regularize  - [1|0] -> [yes|no]
%   lambda      - regularization constant for weight decay. Makes
%		  logistic regression into a support vector machine
%		  for large lambda (cf. Clay Spence). Defaults to 10^-6.
%   lambdasearch- [1|0] -> search for optimal regularization constant lambda
%   eigvalratio - if the data does not fill D-dimensional space,
%                   i.e. rank(x)<D, you should specify a minimum 
%		    eigenvalue ratio relative to the largest eigenvalue
%		    of the SVD.  All dimensions with smaller eigenvalues
%		    will be eliminated prior to the discrimination. 
%


[N,D]=size(x);
iter=0; showcount=0;

if nargin<3 | isempty(v), v=zeros(D,1); vth=0; else vth=v(D+1); v=v(1:D); end;
if nargin<4 | isempty(show); show=0; end;
if nargin<5 | isempty(regularize); regularize=0; end
if nargin<6 | isempty(lambda); lambda=eps; end;
if nargin<7 | isempty(lambdasearch), lambdasearch=0; end
if nargin<8 | isempty(eigvalratio); eigvalratio=0; end;
if regularize, lambda=1e-6; end

% subspace reduction if requested - electrodes might be linked
if eigvalratio
  [U,S,V] = svd(x,0);                        % subspace analysis
  V = V(:,find(diag(S)/S(1,1)>eigvalratio)); % keep significant subspace
  x = x*V;       % map the data to that subspace
  v = V'*v;      % reduce initialization to the subspace
  [N,D]=size(x); % less dimensions now
end

% combine threshold computation with weight vector.
x = [x ones(N,1)];
v = [v; vth];

vold=ones(size(v));

if regularize, lambda = [0.5*lambda*ones(1,D) 0]'; end

% clear warning as we will use it to catch conditioning problems
lastwarn('');

% If lambda increases, the maximum number of iterations is increased from
% maxiter to maxiterlambda
maxiter=100; maxiterlambda=1000;
singularwarning=0; lambdawarning=0; % Initialize warning flags

% IRLS for binary classification of experts (bernoulli distr.)
while ((subspace(vold,v)>1e-7)&(iter<=maxiter)&(~singularwarning)&(~lambdawarning))|iter==0,
    
    
    vold=v;
    mu = bernoull(1,x*v);   % recompute weights
    w = mu.*(1-mu); 
    e = (y - mu);
    grad = x'*e; % - lambda .* v;
    if regularize,
      % lambda=(v'*v)./2;
      % grad=grad-lambda;
      
      grad=grad - lambda .* v; 
    end
    %inc = inv(x'*diag(w)*x+eps*eye(D+1)) * grad;
    inc = inv(x'*(repmat(w,1,D+1).*x)+diag(lambda)*eye(D+1)) * grad;
    
    if strncmp(lastwarn,'Matrix is close to singular or badly scaled.',44)
        warning('Bad conditioning. Suggest reducing subspace.')
        singularwarning=1;
    end
    
    if (norm(inc)>=1000)&regularize, 
        if ~lambdasearch, warning('Data may be perfectly separable. Suggest increasing regularization constant lambda'); end
        lambdawarning=1; 
    end; 
    
    
    % avoid funny outliers that happen with inv    
    if (norm(inc)>=1000)&regularize&lambdasearch, 
        % Increase regularization constant lambda
        lambda=sign(lambda).*abs(lambda.^(1/1.02));
        lambdawarning=0;
        
        if printoutput,
            fprintf('Bad conditioning.  Data may be perfectly separable.  Increasing lambda to: %6.2f\n',lambda(1));
        end
        
        maxiter=maxiterlambda;
        
    elseif (~singularwarning)&(~lambdawarning), 
        % update
        v = v + inc; 
        if sum(isnan(v))
          LL = NaN;
          v = ones(1,size(x,2))*NaN;
          return;
        end
        % exit if converged
        if subspace(vold,v)<1e-7,
            
            
        end;
    end
    
    % exit if taking too long
    
    iter=iter+1;
    if iter>maxiter,
      
    end;
end;

% log-likelihood of the data
LL = sum( x*v.*y) - sum( log(1+ exp(x*v))); 

if eigvalratio
  v = [V*v(1:D);v(D+1)]; % the result should be in the original space
end




return;



function [p]=bernoull(x,eta);
% bernoull() - Computes Bernoulli distribution of x for "natural parameter" eta.
%
% Usage:
%   >> [p] = bernoull(x,eta)
%
% The mean m of a Bernoulli distributions relates to eta as,
% m = exp(eta)/(1+exp(eta));

p = exp(eta.*x - log(1+exp(eta)));

% p = bernoull(1,[x 1]*v);

return;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EXAMPLE OF CLASSIFICATION WITH LOGIT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 
logitp = @(b,x) exp(b(end)+x*b(1:end-1))./(1+exp(b(end)+x*b(1:end-1)));
% p(x;w) = logit
% log( p/(1-p) ) is linear in x
% the decision boundary is x*b = 0
% the class probability depends on the distance from the boundary, whchi is b(0)/norm(b) + x*b/norm(b)
% 

M = [2,4,1]';  % model parameters
X = randn(1000,2);
% model data
Y = logitp(M,X);
r = randn(size(Y))*0.1;
Y = Y+r;
Y2 = zeros(size(Y));
Y2(find(Y>0.5)) = 1;
[class,o] = sort(Y2);
figure(1);clf;
plot(Y(o))

% logist function adds  offset by default, hence we enter only the variable params
[W,~,LL] = ck_decoding_logist(X,Y2,[],0,1,0.01,0,0)
Ypred = X*W(1:2)+W(3);
[~,~,Az] = ck_decoding_auc(Y2,Ypred)

hold on
plot(Ypred(o),'r');

X3 = ones(1000,3);
X3(:,[1,2]) = X;
% log L of data under true model
LLdata = sum( X3*M.*Y2) - sum( log(1+ exp(X3*M)))



