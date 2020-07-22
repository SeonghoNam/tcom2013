function Rate=Proposed(H,N_MaxSupUser,SNR)


dim=size(H);
N_Rx=dim(1); N_Tx=dim(2); N_TotUser=dim(3);

TotUserSet=1:N_TotUser;
InitRateSet=zeros(N_TotUser,1);
PrvSingularValue=zeros(N_TotUser,1);

for iUser=1:N_TotUser
    [u s v]=svd(H(:,:,iUser));
    
    InitRateSet(iUser,iUser)=log2(1+1/N_Rx*det(H(:,:,iUser)*H(:,:,iUser)'))./AvgRateSet(iUser);
%     InitRateSet(:,iUser)=BDRate(H,iUser,1)./AvgRateSet;
    RowSpace(:,:,iUser)=v(:,1:N_Rx);
    
end

[Rate InitUser]=max(sum(InitRateSet,1));
Z{InitUser}=eye(N_Tx);
PrvSingularValue(InitUser)=det(H(:,:,InitUser)*H(:,:,InitUser)');
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
    
    
    ProdEstEig=zeros(N_TotUser,1);
    for iUser=RemaingUserSet


        HeffNew=H(:,:,iUser)*NullspaceBasis;
        ProdEstEig(iUser,iUser)=log2((1/(N_Rx*iSupUser))^N_Rx*(det(HeffNew*HeffNew')));
          
%         M=(RowSpace(:,:,iUser)'*NullspaceBasis);

           
%         cosAng(iUser)=(det(M*M'));
        
        for jUser=SelectedUserSet   
            [u s v]=svd(H(:,:,jUser)*Z{jUser});
            B=v(:,1:N_Rx);
            G=null(H(:,:,iUser)*Z{jUser});
            M=B'*G;
            cosAng=(det(M*M'));
            
            ProdEstEig(jUser,iUser)=PrvSingularValue(jUser)*cosAng;
            SINR=(1/(N_Rx*iSupUser))^N_Rx*ProdEstEig(jUser,iUser);
            ProdEstEig(jUser,iUser)=log2(SINR);
            
%             Heff=H(:,:,jUser)*Z{jUser}*G;
%             [PrvSingularValue(jUser)*cosAng det(Heff*Heff')]
%             ProdEstEig(iUser)=ProdEstEig(iUser)*PrvSingularValue(jUser)*cosAng;
        end

        ProdEstEig(:,iUser)=ProdEstEig(:,iUser)./AvgRateSet;
    end

    
    [tempValue CandidateUser]=max(sum(ProdEstEig,1));
    
%     if log2(((SNR/((iSupUser)*N_Rx))^N_Rx )*tempValue) < log2(((SNR/((iSupUser-1)*N_Rx))^N_Rx )*PrvMetric)
%         break;
%     end

  
    % Update Precoder
    Z{CandidateUser}=NullspaceBasis;
    PrvSingularValue(CandidateUser)=(det(H(:,:,CandidateUser)*Z{CandidateUser}*(H(:,:,CandidateUser)*Z{CandidateUser})'));
    for jUser=SelectedUserSet
        G=null(H(:,:,CandidateUser)*Z{jUser});
        Z{jUser}=Z{jUser}*G;
        PrvSingularValue(jUser)=det((H(:,:,jUser)*Z{jUser})*(H(:,:,jUser)*Z{jUser})');
        

    end
    
    SelectedUserSet=[SelectedUserSet CandidateUser];
    
   
end

% UserSet=SelectedUserSet;
Rate=zeros(N_TotUser,1);
Rate(SelectedUserSet)=c_based(H(:,:,SelectedUserSet),N_MaxSupUser,AvgRateSet(SelectedUserSet));


