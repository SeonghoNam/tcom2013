function Rate=c_based(H,N_MaxSupUser,SNR)


dim=size(H);
N_Rx=dim(1); N_Tx=dim(2); N_TotUser=dim(3);

TotUserSet=1:N_TotUser;
RemaingUserSet=TotUserSet;
for iUser=1:N_TotUser
    InitRate(iUser)=BDRate(H,iUser,SNR);
end

[Rate InitUser]=max(InitRate);
SelectedUserSet=InitUser;

for iSupUser=2:N_MaxSupUser
    RemaingUserSet=TotUserSet;
    RemaingUserSet(SelectedUserSet)=[];
    
    tempRate2=zeros(N_TotUser,1);
    for iUser=RemaingUserSet
        tempRate2(iUser)=BDRate(H,[SelectedUserSet iUser],SNR);
    end
    [tempRate CandidateUser]=max(tempRate2);
    
    if(tempRate>Rate)
        Rate=tempRate;
        SelectedUserSet=[SelectedUserSet CandidateUser];
    else
        break;
    end
end


