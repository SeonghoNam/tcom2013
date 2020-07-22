function RateSet=ER_based(H,N_MaxSupUser,AvgRateSet)


dim=size(H);
N_Rx=dim(1); N_Tx=dim(2); N_TotUser=dim(3);

TotUserSet=1:N_TotUser;

H_tilde=zeros(N_Rx,N_Tx,N_TotUser);

for iUser=1:N_TotUser
    InitRateSet(:,iUser)=BDRate(H,iUser,1)./AvgRateSet;

    [tt b]=rq(H(:,:,iUser));
    H_tilde(:,:,iUser)=b(1:N_Rx,:);
end

[aa InitUser]=max(sum(InitRateSet,1));
[RateSet PrvEigV]=BDRate(H,InitUser,1);
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
    

    tempRateSet=zeros(N_TotUser,N_TotUser);

    for iUser=RemaingUserSet

%         cosAng2=prod(sqrt(diag(NullspaceBasis'*H_tilde(:,:,iUser)'*H_tilde(:,:,iUser)*NullspaceBasis)));
        Heff=H(:,:,iUser)*NullspaceBasis;
        

        M=(NullspaceBasis'*H_tilde(:,:,iUser)')';
        cosAng2=det(M'*M);
        
        EstimatedEigV=PrvEigV*(cosAng2);
        
        NewEigV=svd(Heff);
        eigV=[EstimatedEigV; NewEigV];
        
%         [a realEigV]=BDRate(H,[SelectedUserSet iUser],1);
%         [realEigV eigV [PrvEigV; NewEigV]]
%         realEigV./[PrvEigV; NewEigV]
        
        PowerAllo = WaterFilling_alg(1,eigV,1);
%         PowerAllo = 1/length(eigV)*ones(1,length(eigV));
        

        tempRate=(log2(1+(eigV.^2).*PowerAllo'));
        
        for iNuser=1:iSupUser
            if iNuser==iSupUser
                jUser=iUser;
            else
                jUser=SelectedUserSet(iNuser);
            end
            tempRateSet(jUser,iUser)=sum(tempRate(N_Rx*(iNuser-1)+1:N_Rx*(iNuser)))/AvgRateSet(jUser);
        end

        
        
    end

    
    [tempRate CandidateUser]=max(sum(tempRateSet,1));

    
    [tempRateSet PrvEigV]=BDRate(H,[SelectedUserSet CandidateUser],1);
    
    if(sum(tempRateSet,1)./AvgRateSet>sum(RateSet,1)./AvgRateSet)
        SelectedUserSet=[SelectedUserSet CandidateUser];
        RateSet=tempRateSet;
    else
        break;
    end
end

RateSet=BDRate(H,SelectedUserSet,1);


