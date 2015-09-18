k = 10;
B = basisFactory.makeNonlinearRaisedCos(2,1,[1 5],2);
plot(B.B)
X = Mech.all;
C = any(~isnan(Mech.all),2);
X(isnan(X))=0;
neuralWord(~C)= 0;


[Xout,dm] = buildDesignMatrix(X,neuralWord,'bsStim',B,'deriv',1);
idx = crossvalind('Kfold',length(neuralWord),k);
yhat = nan(size(neuralWord));
t = struct;
parfor ii = 1:k
    yhat = nan(size(neuralWord));
    [b,dev,stats] = glmfit(Xout(idx~=ii,:),neuralWord(idx~=ii),'binomial');
    yhat(idx ==ii) = glmval(b,Xout(idx==ii,:),'logit');
    t(ii).y = yhat;
end
yOut = [];
for ii=1:k
    ydum= t(ii).y;
    ydum(isnan(ydum))=0;
    yOut = [yOut ydum];
end
yOut = sum(yOut,2);

netOut = cvglmnet(Xout,neuralWord,'binomial',[],[],10,[],1);