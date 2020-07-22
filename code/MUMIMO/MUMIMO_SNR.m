clc; clear all; close all;

addpath('.\Gram-Schmidt Process');
N_Tx=8;
N_Rx=2;

N_Algorithm=5;

N_TotUser=20;
N_MaxSupUser=ceil(N_Tx/N_Rx);

TotUserSet=1:N_TotUser;

SNR_dB=0:5:20;
SNR_real=10.^(SNR_dB*0.1);

N_Iteration=3000;

% for iSupUser=1:N_MaxSupUser
%     TotalCombination{iSupUser}=pick(TotUserSet,iSupUser,'');
% end

SumRate=zeros(length(SNR_dB),N_Algorithm);

for iSNR=1:length(SNR_dB)
    SNR_dB(iSNR)
    SNR=SNR_real(iSNR);
    
    for iIter=1:N_Iteration

        H=sqrt(0.5)*(randn([N_Rx,N_Tx,N_TotUser])+1j*randn([N_Rx,N_Tx,N_TotUser]));

        for iAlg=1:N_Algorithm
            switch iAlg
                case 1
%                     disp('Optimal');
%                     tic
                    SumRate(iSNR,iAlg)=SumRate(iSNR,iAlg)+Proposed(H,N_MaxSupUser,SNR);
%                     toc
                case 2
%                     disp('c_based');
%                     tic
                    SumRate(iSNR,iAlg)=SumRate(iSNR,iAlg)+c_based(H,N_MaxSupUser,SNR);
%                     toc
                case 3
%                     disp('prop');
%                     tic
                    SumRate(iSNR,iAlg)=SumRate(iSNR,iAlg)+SRN_based(H,N_MaxSupUser,SNR);
%                     toc
                case 4
%                     disp('ChodalDistance');
%                     tic
                    SumRate(iSNR,iAlg)=SumRate(iSNR,iAlg)+cd_based(H,N_MaxSupUser,SNR);
%                     toc
                case 5
%                     disp('n_based');
%                     tic
                    SumRate(iSNR,iAlg)=SumRate(iSNR,iAlg)+Proposed(H,N_MaxSupUser,SNR);
%                     toc
            end
        end

    end
end

SumRate=SumRate/N_Iteration
% SumRate=real(SumRate);

plot(SNR_dB,SumRate,'-o');
legend('Proposed','c-algorithm','SRN-algorithm','cd-algorithm','Proposed (simple)');
xlabel('SNR [dB]'); ylabel('Sum Capacity [bps/Hz]'); 
