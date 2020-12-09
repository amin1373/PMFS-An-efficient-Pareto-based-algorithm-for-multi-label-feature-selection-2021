function w = mRMMR(BiObj,nSort)
%minimize redundency maximize mean relevance

%Non-dominated sorting
[FrontNo,MaxFNo]=NDSort(BiObj,nSort);

%Calculate the crowding distance of each solution
d=CrowdingDistance(BiObj,FrontNo);
w = FrontNo' + 1./(1+d);

end

