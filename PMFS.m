
function Feature_Ranks= PMFS(X_train ,Y_train)

% X_train = features matrix
% Y_train = labels matrix

[~, feature_num] = size( X_train );
[~, num_class] = size( Y_train );
Lambda = 10;
m1 = pinv(X_train'*X_train + Lambda*eye(feature_num)) * X_train'*Y_train;
m11=max(max(m1,[],2))-m1;
m2=pdist2(m1,m1);
f1=min(m11,[],2);
f2=-sum(m2,2);
BiObj=[f1,f2];
w = mRMMR(BiObj,feature_num);
[v rank]=sort(w);

end