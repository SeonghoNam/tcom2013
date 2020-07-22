function UserSet=Proposed_LimitedFeedback(CQI,CDI,N_MaxSupUser,SNR)


dim=size(CDI);
N_Rx=dim(1); N_Tx=dim(2); N_TotUser=dim(3);

TotUserSet=1:N_TotUser;
PrvSingularValue=zeros(N_TotUser,1);

[tProd InitUser]=max(CQI);
Z{InitUser}=eye(N_Tx);
PrvSingularValue(InitUser)=tProd;
SelectedUserSet=InitUser;

for iSupUser=2:N_MaxSupUser
    RemaingUserSet=TotUserSet;
    RemaingUserSet(SelectedUserSet)=[];
    
   
    if (iSupUser-1)==1
        AggH=CDI(:,:,SelectedUserSet);
    else
        AggH=reshape(permute(CDI(:,:,SelectedUserSet),[2 1 3]),[N_Tx N_Rx*(iSupUser-1)]).';
    end
    
    NullspaceBasis=null(AggH);
    
    
    ProdEstEig=zeros(N_TotUser,1);
    for iUser=RemaingUserSet


        HeffNew=CDI(:,:,iUser)*NullspaceBasis;
        ProdEstEig(iUser)=CQI(iUser)*(det(HeffNew*HeffNew'));
        
%         M=(RowSpace(:,:,iUser)'*NullspaceBasis);

           
%         cosAng(iUser)=(det(M*M'));
        
        for jUser=SelectedUserSet   
            [u s v]=svd(CDI(:,:,jUser)*Z{jUser});
            B=v(:,1:N_Rx);
            G=null(CDI(:,:,iUser)*Z{jUser});
            M=B'*G;
            cosAng(iUser,jUser)=(det(M*M'));
            ProdEstEig(iUser)=ProdEstEig(iUser)*PrvSingularValue(jUser)*cosAng(iUser,jUser);
        end


    end

    
    [tempValue CandidateUser]=max(ProdEstEig);
    

    % Update Precoder
    Z{CandidateUser}=NullspaceBasis;
    PrvSingularValue(CandidateUser)=CQI(CandidateUser)*(det(CDI(:,:,CandidateUser)*Z{CandidateUser}*(CDI(:,:,CandidateUser)*Z{CandidateUser})'));
    for jUser=SelectedUserSet
        G=null(CDI(:,:,CandidateUser)*Z{jUser});
%         [u s v]=svd(CDI(:,:,CandidateUser)*Z{jUser});
%         B=v(:,1:N_rx);
%         M=B'*G;
%         cosAng=(det(M*M'));
        Z{jUser}=Z{jUser}*G;
        PrvSingularValue(jUser)=PrvSingularValue(jUser)*cosAng(CandidateUser,jUser);
        

    end
    
    SelectedUserSet=[SelectedUserSet CandidateUser];
    
   
end

% Rate=log2(((SNR/(length(SelectedUserSet)*N_Rx))^N_Rx )*PrvMetric);
% for jUser=SelectedUserSet
%     Rate=Rate+log2( ((SNR/(length(SelectedUserSet)*N_Rx))^N_Rx )*PrvSingularValue(jUser));
% end

UserSet=SelectedUserSet;
% Rate=c_based(H(:,:,SelecteedUserSet),N_MaxSupUser,SNR);


