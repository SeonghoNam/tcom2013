function UserSet=CD_based(H,N_MaxSupUser,SNR)


dim=size(H);
N_Rx=dim(1); N_Tx=dim(2); N_TotUser=dim(3);

TotUserSet=1:N_TotUser;
for iUser=1:N_TotUser
    InitNorm(iUser)=norm(H(:,:,iUser),'fro')^2;
    

    [tt b]=rq(H(:,:,iUser));
    H_tilde(:,:,iUser)=b(1:N_Rx,:);

end

[tnorm InitUser]=max(InitNorm);
Rate=BDRate(H,InitUser,SNR);

U_tilde=H_tilde(:,:,InitUser);


SelectedUserSet=InitUser;

for iSupUser=2:N_MaxSupUser
    RemaingUserSet=TotUserSet;
    RemaingUserSet(SelectedUserSet)=[];
    
   
    ChordalD=zeros(N_TotUser,1);
    for iUser=RemaingUserSet
        ChordalD(iUser)=N_Rx-trace(H_tilde(:,:,iUser)*U_tilde'*U_tilde*H_tilde(:,:,iUser)');              
    end
   
    [tempCD CandidateUser]=max(ChordalD);
    
    SelectedUserSet=[SelectedUserSet CandidateUser];
    A=H(:,:,CandidateUser)-H(:,:,CandidateUser)*U_tilde'*U_tilde;
    [tt b]=rq(A);
    U_tilde=[U_tilde' b(1:N_Rx,:)']';

%     tempRate=BDRate(H,[SelectedUserSet CandidateUser],SNR);
    
%     if(tempRate>Rate)
%         SelectedUserSet=[SelectedUserSet CandidateUser];
%         Rate=tempRate;
%         
%         A=H(:,:,CandidateUser)-H(:,:,CandidateUser)*U_tilde'*U_tilde;
%         [tt b]=rq(A);
%         U_tilde=[U_tilde' b(1:N_Rx,:)']';
%         
%         
%     else
%         break;
%     end
end
UserSet=SelectedUserSet;
% Rate=c_based(H(:,:,SelectedUserSet),N_MaxSupUser,SNR);

