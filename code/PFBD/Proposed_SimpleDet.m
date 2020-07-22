function Rate=Proposed_SimpleDet(H,N_MaxSupUser,AvgRateSet)


dim=size(H);
N_Rx=dim(1); N_Tx=dim(2); N_TotUser=dim(3);

TotUserSet=1:N_TotUser;
InitRateSet=zeros(N_TotUser,N_TotUser);

for iUser=1:N_TotUser
    [u s v]=svd(H(:,:,iUser));
    
    InitRateSet(iUser,iUser)=log2(1+1/N_Rx*det(H(:,:,iUser)*H(:,:,iUser)'))./AvgRateSet(iUser);
%     InitRateSet(:,iUser)=BDRate(H,iUser,1)./AvgRateSet;
    RowSpace(:,:,iUser)=v(:,1:N_Rx);
    
end

[Rate InitUser]=max(sum(InitRateSet,1));
Z{InitUser}=eye(N_Tx);
PrvSingularValue(InitUser)=det(H(:,:,InitUser)*H(:,:,InitUser)');

% [u s v]=svd(H(:,:,InitUser));
% PrvPrecoderBasis{InitUser}=(v(:,1:N_Rx));

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
    
    
    ProdEstEig=zeros(N_TotUser,N_TotUser);
    for iUser=RemaingUserSet


        HeffNew=H(:,:,iUser)*NullspaceBasis;
        ProdEstEig(iUser,iUser)=log2((1/(N_Rx*iSupUser))^N_Rx*(det(HeffNew*HeffNew')));
        
        M=(RowSpace(:,:,iUser)'*NullspaceBasis);   
        cosAng=(det(M*M'));
        
        for jUser=SelectedUserSet         
            ProdEstEig(jUser,iUser)=PrvSingularValue(jUser)*cosAng;
            SINR=(1/(N_Rx*iSupUser))^N_Rx*ProdEstEig(jUser,iUser);
            ProdEstEig(jUser,iUser)=log2(SINR);
            
%             if(SINR>1)
%                 ProdEstEig(jUser,iUser)=log2(SINR);
%             else
%                 ProdEstEig(jUser,iUser)=SINR*log2(exp(1));
%             end
  

        end

        
        ProdEstEig(:,iUser)=ProdEstEig(:,iUser)./AvgRateSet;

    end

    
    [tempValue CandidateUser]=max(sum(ProdEstEig,1));

    % Update Precoder
    Z{CandidateUser}=NullspaceBasis;
    PrvSingularValue(CandidateUser)=(det(HeffNew*HeffNew'));
    for jUser=SelectedUserSet
        G=null(H(:,:,CandidateUser)*Z{jUser});
        Z{jUser}=Z{jUser}*G;
        PrvSingularValue(jUser)=det((H(:,:,jUser)*Z{jUser})*(H(:,:,jUser)*Z{jUser})');

    end
    
    SelectedUserSet=[SelectedUserSet CandidateUser];
    
   
end

Rate=zeros(N_TotUser,1);
Rate(SelectedUserSet)=c_based(H(:,:,SelectedUserSet),N_MaxSupUser,AvgRateSet(SelectedUserSet));


