function [pred_r, null_pred_r, pred_x, pred_y, brainWei] = func_cca_pred_analysis_LOOCV(X, Y)

num_perm = 1000;
num_dim = min(size(X,2), size(Y,2));
num_obsv = size(X, 1);

null_pred_r = [];
pred_r = [];

pred_x=[];
pred_y=[];

for i=1:num_obsv
    Xtrain = X;
    Xtest = X(i,:);
    Xtrain(i,:) = [];
    
    Ytrain = Y;
    Ytest = Y(i,:);
    Ytrain(i,:) = [];
    
    [A,B]=canoncorr(Xtrain, Ytrain);
    pred_x(i,:) = Xtest*A;
    brainWei(:,i)=A(:,1);
    pred_y(i,:) = Ytest*B;
end

for i = 1:size(pred_x,2)
  pred_r(i) = corr(pred_x(:,i), pred_y(:,i), 'type', 'Pearson');
end

for inull = 1:num_perm
  null_pred_x = [];
  null_pred_y = [];

  null_idx = randperm(num_obsv);
  null_X = X(null_idx, :);
  for i = 1:num_obsv
	Xtrain = null_X;
    Xtest = null_X(i,:);
    Xtrain(i,:) = [];
    
    Ytrain = Y;
    Ytest = Y(i,:);
    Ytrain(i,:) = [];
    
    [A,B] = canoncorr(Xtrain, Ytrain);
    
    null_pred_x(i,:) = Xtest*A;
    null_pred_y(i,:) = Ytest*B;
  end
  for i = 1:size(null_pred_x,2)
    null_pred_r(inull, i) = corr(null_pred_x(:,i), null_pred_y(:,i), 'type', 'Pearson');
  end
end

