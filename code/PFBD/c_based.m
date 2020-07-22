function RateSet=c_based(H,N_MaxSupUser,AvgRateSet)


dim=size(H);
N_Rx=dim(1); N_Tx=dim(2); N_TotUser=dim(3);

TotUserSet=1:N_TotUser;
for iUser=1:N_TotUser
    InitRateSet(:,iUser)=BDRate(H,iUser,1)./AvgRateSet;
end

[Rate InitUser]=max(sum(InitRateSet,1));
RateSet=InitRateSet(:,InitUser);
SelectedUserSet=InitUser;

for iSupUser=2:N_MaxSupUser
    RemaingUserSet=TotUserSet;
    RemaingUserSet(SelectedUserSet)=[];
    
    tempRateSet=zeros(N_TotUser,N_TotUser);
    for iUser=RemaingUserSet
        tempRateSet(:,iUser)=BDRate(H,[SelectedUserSet iUser],1)./AvgRateSet;
    end
    [tempRate CandidateUser]=max(sum(tempRateSet,1));
    
    
    if(tempRate>Rate)
        Rate=tempRate;
        RateSet=tempRateSet(:,CandidateUser);
        SelectedUserSet=[SelectedUserSet CandidateUser];
    else
        break;
    end
end

RateSet=BDRate(H,SelectedUserSet,1);
