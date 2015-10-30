k = 10;
doNormal =0;
noProx = 1;
B = basisFactory.makeNonlinearRaisedCos(4,1,[0 6],1);
B2 = basisFactory.makeNonlinearRaisedCos(3,1,[1 10],5);

plot(B2.B)
plot(B.B)
drawnow

X = Mech.all;
X = [Mech.Fx sqrt(Mech.Fy.^2+Mech.Fz.^2) unwrap(atan2(Mech.Fz,Mech.Fy)) Mech.Mx sqrt(Mech.My.^2+Mech.Mz.^2) unwrap(atan2(Mech.Mz,Mech.My))];
C = any(~isnan(X),2);
X(isnan(X))=0;
neuralWord(~C)= 0;

if noProx
    X = X(~prox,:);
    neuralWord = neuralWord(~prox);
    C = C(~prox);
end

[Xout,dm] = buildDesignMatrix(X,neuralWord,'bsStim',B,'deriv',0,'hist',1);%,'bsSpike',B2);
% Xout = [Xout [0;neuralWord(1:end-1)]];
idx = crossvalind('Kfold',length(neuralWord),k);
idxTest = idx(C);
t = struct;
coeffs = struct;

% [b,dev,stats] = glmfit(Xout,neuralWord,'binomial');

% netOut = cvglmnet(Xout,neuralWord,'binomial',[],[],10,[],1);




parfor ii = 1:k
    yhat = nan(size(neuralWord(find(C==1))));
%     yhat = nan(size(neuralWord));
    yhat_normal = yhat;
    [b,dev,stats] = glmfit(Xout(idx~=ii,:),neuralWord(idx~=ii),'binomial');
    Xtest = Xout(find(C),:);
    yhat(idxTest==ii) = glmval(b,Xtest(idxTest==ii,:),'logit');
    
    t(ii).y = yhat;
    t(ii).b = b;
    %     coeffs(ii).coef = b;
    %     if doNormal
    %         yhat_normal(idxTest==ii) = glmval(b([1 B2.edim+2:end]),Xtest(idxTest==ii,[B2.edim+1:end]),'identity');
    %         t(ii).yNormal = yhat_normal;
    %
    %     end
end
% [XoutHist,dm] = buildDesignMatrix(X,neuralWord,'bsStim',B,'deriv',0,'hist',1,'bsSpike',B2);
% 
% b2 = glmfit(XoutHist,neuralWord,'binomial');
% yhatNormal = glmval(b2([1 B2.edim+2:end]),XoutHist(:,[B2.edim+1:end]),'identity');
% spikeWeights = buildGLM.combineWeights(dm,b2(2:end));
% spikeWeights = spikeWeights.hist;
% fprintf('Simulating')
% [simOut,lambda] = simGLM4(yhatNormal,spikeWeights.data,500);
% 
% yOut = [];
% if doNormal
%     yOut_normal = [];
% end
% for ii=1:k
%     ydum= t(ii).y;
%     ydum(isnan(ydum))=0;
%     yOut = [yOut ydum];
%     if doNormal
%         ydum = t(ii).yNormal;
%         ydum(isnan(ydum))=0;
%         yOut_normal = [yOut_normal ydum];
%     end
%     
% end
% if doNormal
%     yOut_normal = sum(yOut_normal,2);
% end
% yOut = sum(yOut,2);
% yTest = neuralWord(C);
