function Rate=SRN_based(H,N_MaxSupUser,SNR)


dim=size(H);
N_Rx=dim(1); N_Tx=dim(2); N_TotUser=dim(3);

TotUserSet=1:N_TotUser;
InitNorm=zeros(N_TotUser,1);

for iUser=1:N_TotUser
    InitNorm(iUser)=prod(diag(H(:,:,iUser)*H(:,:,iUser)'));
    

end

[tnorm InitUser]=max(InitNorm);
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
    
    
    SquaredRowNorm=zeros(N_TotUser,1);
    for iUser=RemaingUserSet


        HeffNew=H(:,:,iUser)*NullspaceBasis;
        SquaredRowNorm(iUser)=prod(diag(HeffNew*HeffNew'));
        for jUser=SelectedUserSet
            G=null(H(:,:,iUser)*Z{jUser});
            HeffPrv=H(:,:,jUser)*Z{jUser}*G;
            SquaredRowNorm(iUser)=SquaredRowNorm(iUser)*prod(diag(HeffPrv*HeffPrv'));
        end


    end

    
    [tempValue CandidateUser]=max(SquaredRowNorm);
    
    % Update Precoder
    Z{CandidateUser}=NullspaceBasis;
    for jUser=SelectedUserSet
        G=null(H(:,:,CandidateUser)*Z{jUser});
        Z{jUser}=Z{jUser}*G;
    end
    
    SelectedUserSet=[SelectedUserSet CandidateUser];
    
   
end

Rate=c_based(H(:,:,SelectedUserSet),N_MaxSupUser,SNR);


