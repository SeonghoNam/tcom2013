function Rate=SRN_based(H,N_MaxSupUser,AvgRateSet)


dim=size(H);
N_Rx=dim(1); N_Tx=dim(2); N_TotUser=dim(3);

TotUserSet=1:N_TotUser;
InitRateSet=zeros(N_TotUser,N_TotUser);

for iUser=1:N_TotUser
    InitRateSet(iUser,iUser)=sum(log2(1+1/N_Rx*(diag(H(:,:,iUser)*H(:,:,iUser)'))))./AvgRateSet(iUser);
end

[aa InitUser]=max(sum(InitRateSet,1));

Z{InitUser}=eye(N_Tx);

SelectedUserSet=InitUser;

for iSupUser=2:N_MaxSupUser
    RemaingUserSet=TotUserSet;
    RemaingUserSet(SelectedUserSet)=[];
    
   
    if (iSupUser-1)==1
        AggH=H(:,:,SelectedUserSet);
    else
        AggH=reshape(permute(H(:,:,SelectedUserSet),[2 1 3]),[N_Tx N_Rx*(iSupUser-1)]).';
    end
    
    NullspaceBasis=null(AggH);
    
    
    SquaredRowNorm=zeros(N_TotUser,N_TotUser);
    for iUser=RemaingUserSet


        HeffNew=H(:,:,iUser)*NullspaceBasis;
        SquaredRowNorm(iUser,iUser)=sum(log2(1+1/(N_Rx*iSupUser)*diag(HeffNew*HeffNew')));
        for jUser=SelectedUserSet
            G=null(H(:,:,iUser)*Z{jUser});
            HeffPrv=H(:,:,jUser)*Z{jUser}*G;
            SquaredRowNorm(jUser,iUser)=sum(log2(1+1/(N_Rx*iSupUser)*diag(HeffPrv*HeffPrv')));
        end

        SquaredRowNorm(:,iUser)=SquaredRowNorm(:,iUser)./AvgRateSet;
    end

    
    [tempValue CandidateUser]=max(sum(SquaredRowNorm,1));
    
    % Update Precoder
    Z{CandidateUser}=NullspaceBasis;
    for jUser=SelectedUserSet
        G=null(H(:,:,CandidateUser)*Z{jUser});
        Z{jUser}=Z{jUser}*G;
    end
    
    SelectedUserSet=[SelectedUserSet CandidateUser];
    
   
end

Rate=zeros(N_TotUser,1);
Rate(SelectedUserSet)=c_based(H(:,:,SelectedUserSet),N_MaxSupUser,AvgRateSet(SelectedUserSet));

