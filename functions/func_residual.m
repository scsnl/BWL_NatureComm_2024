function [resid]=func_residual(y, X)

% b = regress(Y, [ones(length(X),1) X]);
% Yhat = [ones(length(X),1) X]*b;
% resid = Y - Yhat;

resid=zeros(size(y));
for d=1:size(y,2)
    [b,bint,resid(:,d)] = regress(y(:,d),[ones(length(X(:,1)),1), X]);
end