function [Rate eigValue]=BDRate(H,UserSet,SNR)

dim=size(H);
N_Rx=dim(1); N_Tx=dim(2); N_ToTUser=dim(3);

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


    
    eigV=[eigV; svd(Heff)];
    

end

PowerAllo = WaterFilling_alg(1,eigV,1/SNR);   

% for iNuser=1:N_SelectedUser
%     aa=N_Rx*(iNuser-1)+1:N_Rx*(iNuser);
%     tempRate(iNuser)=log2( det(eye(length(aa))+SNR*diag(eigV(aa))^2*diag(PowerAllo(aa))) );
% end

% tempRate
tempRate=(log2(1+SNR*(eigV.^2).*PowerAllo'));
RateSet=zeros(N_ToTUser,1);
for iNuser=1:N_SelectedUser
    iUser=UserSet(iNuser);
    RateSet(iUser)=sum(tempRate(N_Rx*(iNuser-1)+1:N_Rx*(iNuser)));
end

% tempRate=log2( det(eye(length(eigV))+SNR*diag(eigV)^2*diag(PowerAllo)) );

if nargout==1
    Rate=RateSet;
elseif nargout==2
    Rate=RateSet;
    eigValue=eigV;
end
    
