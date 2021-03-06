function Rate=Proposed(H,N_MaxSupUser,SNR)


dim=size(H);
N_Rx=dim(1); N_Tx=dim(2); N_TotUser=dim(3);

TotUserSet=1:N_TotUser;
InitProd=zeros(N_TotUser,1);
PrvSingularValue=zeros(N_TotUser,1);

for iUser=1:N_TotUser
    [u s v]=svd(H(:,:,iUser));
    
    InitProd(iUser)=det(H(:,:,iUser)*H(:,:,iUser)');
%     InitProd(iUser)=prod(diag(s))^2;
    RowSpace(:,:,iUser)=v(:,1:N_Rx);
end

[tProd InitUser]=max(InitProd);
Z{InitUser}=eye(N_Tx);
PrvSingularValue(InitUser)=tProd;
PrvMetric=tProd;
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
        ProdEstEig(iUser)=(det(HeffNew*HeffNew'));
        
%         M=(RowSpace(:,:,iUser)'*NullspaceBasis);

           
%         cosAng(iUser)=(det(M*M'));
        
        for jUser=SelectedUserSet   
            [u s v]=svd(H(:,:,jUser)*Z{jUser});
            B=v(:,1:N_Rx);
            G=null(H(:,:,iUser)*Z{jUser});
            M=B'*G;
            cosAng=(det(M*M'));
%             Heff=H(:,:,jUser)*Z{jUser}*G;
%             [PrvSingularValue(jUser)*cosAng det(Heff*Heff')]
            ProdEstEig(iUser)=ProdEstEig(iUser)*PrvSingularValue(jUser)*cosAng;
        end


    end

    
    [tempValue CandidateUser]=max(ProdEstEig);
    
%     if log2(((SNR/((iSupUser)*N_Rx))^N_Rx )*tempValue) < log2(((SNR/((iSupUser-1)*N_Rx))^N_Rx )*PrvMetric)
%         break;
%     end

    PrvMetric=tempValue;
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
Rate=c_based(H(:,:,SelectedUserSet),N_MaxSupUser,SNR);


