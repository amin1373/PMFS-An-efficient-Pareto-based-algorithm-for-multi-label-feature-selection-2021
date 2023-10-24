function [FrontNo,MaxFNo] = NDSort(BiObj, nSort)
%NDSort - Do non-dominated sorting by efficient non-dominated sort.
ObjMatrix=BiObj;
%   FrontNo = NDSort(ObjMatrix, nSort) does non-dominated sorting on ObjMatrix, where ObjMatrix is the
%   matrix of objective values, and nSort is the number
%   of features to be sorted at least. FrontNo(i) denotes the front
%   number of the i-th feature. The features have not been sorted are
%   assigned a front number of inf.
%
%   In particular, nSort = 1 indicates finding only the first non-dominated
%   front, nSort = size(ObjMatrix,1)/2 indicates sorting only half the features
%   , and nSort = inf indicates sorting the whole features.

    [N,M] = size(ObjMatrix);
    if M < 5 || N < 500
        % Use efficient non-dominated sort with sequential search (ENS-SS)
        [FrontNo,MaxFNo] = ENS_SS(ObjMatrix,nSort);
    else
        % Use tree-based efficient non-dominated sort (T-ENS)
        [FrontNo,MaxFNo] = T_ENS(ObjMatrix,nSort);
    end
end

function [FrontNo,MaxFNo] = ENS_SS(PopObj,nSort)
    [PopObj,~,Loc] = unique(PopObj,'rows');
    Table          = hist(Loc,1:max(Loc));
    [N,M]          = size(PopObj);
    [PopObj,rank]  = sortrows(PopObj);
    FrontNo        = inf(1,N);
    MaxFNo         = 0;
    while sum(Table(FrontNo<inf)) < min(nSort,length(Loc))
        MaxFNo = MaxFNo + 1;
        for i = 1 : N
            if FrontNo(i) == inf
                Dominated = false;
                for j = i-1 : -1 : 1
                    if FrontNo(j) == MaxFNo
                        m = 2;
                        while m <= M && PopObj(i,m) >= PopObj(j,m)
                            m = m + 1;
                        end
                        Dominated = m > M;
                        if Dominated || M == 2
                            break;
                        end
                    end
                end
                if ~Dominated
                    FrontNo(i) = MaxFNo;
                end
            end
        end
    end
    FrontNo(rank) = FrontNo;
    FrontNo = FrontNo(:,Loc);
end

function [FrontNo,MaxFNo] = T_ENS(PopObj,nSort)
    [PopObj,~,Loc] = unique(PopObj,'rows');
    Table          = hist(Loc,1:max(Loc));
	[N,M]          = size(PopObj);
    [PopObj,rank]  = sortrows(PopObj);
    FrontNo        = inf(1,N);
    MaxFNo         = 0;
    Forest    = zeros(1,N);
    Children  = zeros(N,M-1);
    LeftChild = zeros(1,N) + M;
    Father    = zeros(1,N);
    Brother   = zeros(1,N) + M;
    [~,ORank] = sort(PopObj(:,2:M),2,'descend');
    ORank     = ORank + 1;
    while sum(Table(FrontNo<inf)) < min(nSort,length(Loc))
        MaxFNo = MaxFNo + 1;
        root   = find(FrontNo==inf,1);
        Forest(MaxFNo) = root;
        FrontNo(root)  = MaxFNo;
        for p = 1 : N
            if FrontNo(p) == inf
                Pruning = zeros(1,N);
                q = Forest(MaxFNo);
                while true
                    m = 1;
                    while m < M && PopObj(p,ORank(q,m)) >= PopObj(q,ORank(q,m))
                        m = m + 1;
                    end
                    if m == M
                        break;
                    else
                        Pruning(q) = m;
                        if LeftChild(q) <= Pruning(q)
                            q = Children(q,LeftChild(q));
                        else
                            while Father(q) && Brother(q) > Pruning(Father(q))
                                q = Father(q);
                            end
                            if Father(q)
                                q = Children(Father(q),Brother(q));
                            else
                                break;
                            end
                        end
                    end
                end
                if m < M
                    FrontNo(p) = MaxFNo;
                    q = Forest(MaxFNo);
                    while Children(q,Pruning(q))
                        q = Children(q,Pruning(q));
                    end
                    Children(q,Pruning(q)) = p;
                    Father(p) = q;
                    if LeftChild(q) > Pruning(q)
                        Brother(p)   = LeftChild(q);
                        LeftChild(q) = Pruning(q);
                    else
                        bro = Children(q,LeftChild(q));
                        while Brother(bro) < Pruning(q)
                            bro = Children(q,Brother(bro));
                        end
                        Brother(p)   = Brother(bro);
                        Brother(bro) = Pruning(q);
                    end
                end
            end
        end
    end
    FrontNo(rank) = FrontNo;
    FrontNo = FrontNo(:,Loc);
end