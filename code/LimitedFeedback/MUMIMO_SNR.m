clc; clear all; close all;

addpath('.\Gram-Schmidt Process');
N_Tx=4;
N_Rx=2;



N_TotUser=20;
N_MaxSupUser=ceil(N_Tx/N_Rx);

TotUserSet=1:N_TotUser;

SNR_dB=0:5:20;
SNR_dB=[0 20];
SNR_real=10.^(SNR_dB*0.1);

N_Iteration=5000;


Qb_CDI=[1 1 1 1];
Qb_CQI=[0 1 2 3];
N_extraAlgorithm=2;
N_Algorithm=length(Qb_CDI)+N_extraAlgorithm;

for iAlg=1:length(Qb_CDI)
    [Level{iAlg} Value{iAlg}]=GetDetQuantizerLevel(1,N_Tx,N_Rx,Qb_CQI(iAlg));
    for iUser=1:N_TotUser
        Codebook{iUser,iAlg}=GenerateCDICodebook(N_Tx,N_Rx,Qb_CDI(iAlg));
    end
end

        
SumRate=zeros(length(SNR_dB),N_Algorithm);
TotalRate=zeros(N_Iteration,N_Algorithm,length(SNR_dB));
for iSNR=1:length(SNR_dB)
    SNR_dB(iSNR)
    SNR=SNR_real(iSNR);
    

%     for iAlg=1:length(Qb_CDI)
%         [Level{iAlg} Value{iAlg}]=GetDetQuantizerLevel(SNR,N_Tx,N_Rx,Qb_CQI(iAlg));
%         for iUser=1:N_TotUser
%             Codebook{iUser,iAlg}=GenerateCDICodebook(N_Tx,N_Rx,Qb_CDI(iAlg));
%         end
%     end

    for iIter=1:N_Iteration

        iIter
        H=sqrt(0.5)*(randn([N_Rx,N_Tx,N_TotUser])+1j*randn([N_Rx,N_Tx,N_TotUser]));
        


        
        for iAlg=1:N_Algorithm
            if iAlg>N_extraAlgorithm
                for iUser=1:N_TotUser
                    CQI(iUser)=DetQuantizer(det(H(:,:,iUser)*H(:,:,iUser)'),Level{iAlg-N_extraAlgorithm},Value{iAlg-N_extraAlgorithm});
%                     CDI(:,:,iUser)=CDIQuantizer(H(:,:,iUser),Codebook{iUser,iAlg-N_extraAlgorithm});
                    DET(iUser)=det(H(:,:,iUser)*H(:,:,iUser)');
                    [a b c]=svd(H(:,:,iUser));
                    CDI(:,:,iUser)=c(:,1:N_Rx)';
% %                     
                end
%                 CQI
            end
            
           if(iAlg==1)
%                continue;
                UserSet=Proposed(H,N_MaxSupUser,SNR);
                TotalRate(iIter,iAlg,iSNR) = BDRate(H,UserSet,SNR);
                SumRate(iSNR,iAlg)=SumRate(iSNR,iAlg)+TotalRate(iIter,iAlg,iSNR);

%                     SumRate(iSNR,iAlg)=SumRate(iSNR,iAlg)+BDRate(H,UserSet,SNR);
           elseif(iAlg==2)
%                continue;
                UserSet=CD_based(H,N_MaxSupUser,SNR);
                TotalRate(iIter,iAlg,iSNR) = BDRate(H,UserSet,SNR);
                SumRate(iSNR,iAlg)=SumRate(iSNR,iAlg)+TotalRate(iIter,iAlg,iSNR);

%                     SumRate(iSNR,iAlg)=SumRate(iSNR,iAlg)+BDRate(H,UserSet,SNR);
           elseif(iAlg==3)
%                 continue;
                UserSet=CD_based(CDI,N_MaxSupUser,SNR);
                TotalRate(iIter,iAlg,iSNR) = BDRate_Feedback(H,CDI,UserSet,SNR);
                SumRate(iSNR,iAlg)=SumRate(iSNR,iAlg)+TotalRate(iIter,iAlg,iSNR);

%                     SumRate(iSNR,iAlg)=SumRate(iSNR,iAlg)+BDRate_Feedback(H,CDI_Real,UserSet,SNR);
% 
           else
                UserSet=Proposed_LimitedFeedback(CQI,CDI,N_MaxSupUser,SNR);
                TotalRate(iIter,iAlg,iSNR) = BDRate_Feedback(H,CDI,UserSet,SNR);
                SumRate(iSNR,iAlg)=SumRate(iSNR,iAlg)+TotalRate(iIter,iAlg,iSNR);

%                     SumRate(iSNR,iAlg)=SumRate(iSNR,iAlg)+BDRate_Feedback(H,CDI,UserSet,SNR);
           end

            
        end

    end
end

SumRate=SumRate/N_Iteration
% SumRate=real(SumRate);

plot(SNR_dB,SumRate,'-o');
legend('Proposed (Perfect CSI)','cd-algorithm (Perfect CSI)','cd-algorithm (Perfect CDI)','Proposed (Perfect CDI with q=0)','Proposed (Perfect CDI with q=1)','Proposed (Perfect CDI with q=2)','Proposed (Perfect CDI with q=3)',2);
xlabel('SNR [dB]'); ylabel('Average Sum Rate [bps/Hz]'); 


for iAlg = 1:N_Algorithm
    [cdf(:,iAlg) x(:,iAlg)] = plotcdf(TotalRate(:,iAlg,1));
end
figure(2)
% cdf(:,1:2)=[]; x(:,1:2)=[]; 
plot(x,cdf); 
ylim([0 1]);
legend('Proposed (Perfect CSI)','cd-algorithm (Perfect CSI)','cd-algorithm (Perfect CDI)','Proposed (Perfect CDI with q=0)','Proposed (Perfect CDI with q=1)','Proposed (Perfect CDI with q=2)','Proposed (Perfect CDI with q=3)',2);

clear cdf, x;
for iAlg = 1:N_Algorithm
    [cdf(:,iAlg) x(:,iAlg)] = plotcdf(TotalRate(:,iAlg,2));
end
figure(3)
% cdf(:,1:2)=[]; x(:,1:2)=[]; 
plot(x,cdf); 
ylim([0 1]);
legend('Proposed (Perfect CSI)','cd-algorithm (Perfect CSI)','cd-algorithm (Perfect CDI)','Proposed (Perfect CDI with q=0)','Proposed (Perfect CDI with q=1)','Proposed (Perfect CDI with q=2)','Proposed (Perfect CDI with q=3)',2);
