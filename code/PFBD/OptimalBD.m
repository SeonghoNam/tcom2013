function Rate=OptimalBD(H,TotalCombination,AvgRateSet)

% dim=size(H);
% N_Rx=dim(1); N_Tx=dim(2); 

for iSupUser=1:length(TotalCombination)
    Comb=TotalCombination{iSupUser};
    for icomb=1:size(Comb,1)
        UserSet=Comb(icomb,:);
        tempRateSet(:,icomb)=BDRate(H,UserSet,1)./AvgRateSet;
    end
    [tempRate2(iSupUser) index]=max(sum(tempRateSet,1));
    tempRateSet2(:,iSupUser)=tempRateSet(:,index);
end

[maxRate index]=max(sum(tempRateSet2,1));
Rate=tempRateSet2(:,index).*AvgRateSet;

% end