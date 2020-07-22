function Rate=ER_based2(H,N_MaxSupUser,AvgRateSet)


dim=size(H);
N_Rx=dim(1); N_Tx=dim(2); N_TotUser=dim(3);

TotUserSet=1:N_TotUser;
InitRateSet=zeros(N_TotUser,N_TotUser);

for iUser=1:N_TotUser
    InitRateSet(:,iUser)=BDRate(H,iUser,1)./AvgRateSet;

end

[Rate InitUser]=max(sum(InitRateSet,1));
Z{InitUser}=eye(N_Tx);
[u s v]=svd(H(:,:,InitUser));
PrvPrecoderBasis{InitUser}=(v(:,1:N_Rx));

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
        ProdEstEig(iUser,iUser)=log2((1/(N_Rx*iSupUser))^N_Rx*prod(eig(HeffNew*HeffNew')));
%         ProdEstEig(iUser,iUser)=sum(log2(1+(1/(N_Rx*iSupUser))*eig(HeffNew*HeffNew')));
        for jUser=SelectedUserSet
            G=null(H(:,:,iUser)*Z{jUser});


            M=(PrvPrecoderBasis{jUser}'*G);
%             cosAng=det(M'*M);
            cosAng=prod(svd(M));
            
%                         
%             
            ProdEstEig(jUser,iUser)=(prod(svd(H(:,:,jUser)*Z{jUser}))*cosAng)^2;
            SINR=(1/(N_Rx*iSupUser))^N_Rx*ProdEstEig(jUser,iUser);
            
            if(SINR>1)
                ProdEstEig(jUser,iUser)=log2(SINR);
            else
                ProdEstEig(jUser,iUser)=SINR*log2(exp(1));
            end
  
%             cosAng=(trace(M'*M)/N_Rx);
%             ProdEstEig(jUser,iUser)=sum(log2(1+(1/(N_Rx*iSupUser))*(svd(H(:,:,jUser)*Z{jUser})).^2*cosAng));
        end

        
        ProdEstEig(:,iUser)=ProdEstEig(:,iUser)./AvgRateSet;

    end

    
    [tempValue CandidateUser]=max(sum(ProdEstEig,1));

    % Update Precoder
    Z{CandidateUser}=NullspaceBasis;
    [u s v]=svd(H(:,:,CandidateUser)*Z{CandidateUser});
    PrvPrecoderBasis{CandidateUser}=(v(:,1:N_Rx));
    for jUser=SelectedUserSet
        G=null(H(:,:,CandidateUser)*Z{jUser});

        Z{jUser}=Z{jUser}*G;
        [u s v]=svd(H(:,:,jUser)*Z{jUser});

        
%         H(:,:,jUser)*Z{jUser}
        PrvPrecoderBasis{jUser}=(v(:,1:N_Rx));

    end
    
    SelectedUserSet=[SelectedUserSet CandidateUser];
    
   
end

Rate=zeros(N_TotUser,1);
Rate(SelectedUserSet)=c_based(H(:,:,SelectedUserSet),N_MaxSupUser,AvgRateSet(SelectedUserSet));


