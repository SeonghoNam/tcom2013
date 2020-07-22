function RateSet=c_based_low(H,N_MaxSupUser,AvgRateSet)
SNR=1;

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

    if (iSupUser-1)==1
        AggH=H(:,:,SelectedUserSet);
    else
        AggH=reshape(permute(H(:,:,SelectedUserSet),[2 1 3]),[N_Tx N_Rx*(iSupUser-1)]).';
    end
    
    
    Precoder=null(AggH);

    tempRate=zeros(N_TotUser,1);
    for iUser=RemaingUserSet
        Heff=H(:,:,iUser)*Precoder;
        eigV=svd(Heff);
        
        PowerAllo = WaterFilling_alg(1/iSupUser,eigV,1/SNR);   
        tempRate(iUser)=log2( det(eye(length(eigV))+SNR*diag(eigV)^2*diag(PowerAllo))) / AvgRateSet(iUser);      

    end
    
    [temp CandidateUser]=max(tempRate);
    
    tempRateSet=BDRate(H,[SelectedUserSet CandidateUser],SNR);
    tempRate=sum(tempRateSet);
    
    if(tempRate>Rate)
        Rate=tempRate;
        RateSet=tempRateSet;
        SelectedUserSet=[SelectedUserSet CandidateUser];
    else
        break;
    end
end

% % SelectedUserSet
% Rate=c_based(H(:,:,SelectedUserSet),N_MaxSupUser,SNR);

% end