function [Rate eigValue]=BDRate(H,UserSet,SNR)

dim=size(H);
N_Rx=dim(1); N_Tx=dim(2); 

N_SelectedUser=length(UserSet);

if N_SelectedUser==1
    AggH=H(:,:,UserSet);
else
    AggH=reshape(permute(H(:,:,UserSet),[2 1 3]),[N_Tx N_Rx*N_SelectedUser]).';
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
    Precoder=null(H_tilde);
    Heff=H(:,:,iUser)*Precoder;

%     [U S V]=svd(Heff);
%     Heff=H(:,:,iUser)*Precoder*V(:,1:rank(Heff));
% 
%     lamda=diag(S); 
    
    eigV=[eigV; svd(Heff)];
    
%     PowerAllo = WaterFilling_alg(1,lamda,1/SNR);   

    
%     tempRate=tempRate+log2( det(eye(N_Rx)+(SNR)*Heff*diag(PowerAllo)*Heff') );
%     tempRate(iNuser)=log2( det(eye(N_Rx)+(SNR/N_Tx)*diag(svd(Heff))^2) );
end

PowerAllo = WaterFilling_alg(1,eigV,1/SNR);   

% for iNuser=1:N_SelectedUser
%     aa=N_Rx*(iNuser-1)+1:N_Rx*(iNuser);
%     tempRate(iNuser)=log2( det(eye(length(aa))+SNR*diag(eigV(aa))^2*diag(PowerAllo(aa))) );
% end

% tempRate
tempRate=sum(log2(1+SNR*(eigV.^2).*PowerAllo'));
% tempRate=log2( det(eye(length(eigV))+SNR*diag(eigV)^2*diag(PowerAllo)) );

if nargout==1
    Rate=tempRate;
elseif nargout==2
    Rate=tempRate;
    eigValue=eigV;
end
    
