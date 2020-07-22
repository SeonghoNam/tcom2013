function [Rate eigValue]=BDRate_Feedback(H,H_hat,UserSet,SNR)

dim=size(H);
N_Rx=dim(1); N_Tx=dim(2); 

N_SelectedUser=length(UserSet);

if N_SelectedUser==1
    AggH=H_hat(:,:,UserSet);
else
    AggH=reshape(permute(H_hat(:,:,UserSet),[2 1 3]),[N_Tx N_Rx*N_SelectedUser]).';
end


% AggH=[];
% for iUser=UserSet
%     AggH=[AggH' H(:,:,iUser)']';
% end



tempRate=0;
eigV=[];
for iNuser=1:N_SelectedUser
    iUser=UserSet(iNuser);
    H_tilde=AggH;
    H_tilde(N_Rx*(iNuser-1)+1:N_Rx*(iNuser),:)=[];
    Precoder(:,:,iUser)=null(H_tilde);
end
for iNuser=1:N_SelectedUser
    iUser=UserSet(iNuser);
    Hint=zeros(N_Rx,N_Rx);
    Hmy=zeros(N_Rx,N_Rx);
    for jNuser=1:N_SelectedUser
        jUser=UserSet(jNuser);
        if iUser==jUser
            Hmy=Hmy+H(:,:,iUser)*Precoder(:,:,jUser);
        else
%             Hmy=Hmy+H(:,:,iUser)*Precoder(:,:,jUser);
            Hint=Hint+H(:,:,iUser)*Precoder(:,:,jUser);
        end
    end

    UserRate(iUser)=log2(real(det(eye(N_Rx)+(SNR/(N_SelectedUser*N_Rx))*Hmy*inv(eye(N_Rx)+(SNR/(N_SelectedUser*N_Rx))*Hint*Hint')*Hmy')));
%     UserRate(iUser)=log2(real(det(eye(N_Rx)+(SNR/(N_SelectedUser*N_Rx))*Hmy*Hmy'))) -log2(real(det(eye(N_Rx)+(SNR/(N_SelectedUser*N_Rx))*Hint*Hint')));
    
end
% UserRate
% sum(UserRate)
tempRate=sum(UserRate);

% PowerAllo = WaterFilling_alg(1,eigV,1/SNR);   
% 
% % for iNuser=1:N_SelectedUser
% %     aa=N_Rx*(iNuser-1)+1:N_Rx*(iNuser);
% %     tempRate(iNuser)=log2( det(eye(length(aa))+SNR*diag(eigV(aa))^2*diag(PowerAllo(aa))) );
% % end
% 
% % tempRate
% tempRate=sum(log2(1+SNR*(eigV.^2).*PowerAllo'));
% tempRate=log2( det(eye(length(eigV))+SNR*diag(eigV)^2*diag(PowerAllo)) );

if nargout==1
    Rate=tempRate;
elseif nargout==2
    Rate=tempRate;
    eigValue=eigV;
end
    
