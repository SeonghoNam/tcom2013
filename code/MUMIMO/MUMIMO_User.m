clc; clear all; close all;

N_Tx=6;
N_Rx=2;

N_Algorithm=4;


N_MaxSupUser=ceil(N_Tx/N_Rx);


SNR_dB=20;
SNR_real=10.^(SNR_dB*0.1);

N_Iteration=1000;



N_User=[5 10 15 20 25 30 35 40]; 
% N_User=[100]; 
SumRate=zeros(length(N_User),N_Algorithm);

for iNuser=1:length(N_User);
    N_TotUser=N_User(iNuser)
    
    TotUserSet=1:N_TotUser;
    
    SNR=SNR_real;

    
    for iIter=1:N_Iteration

        H=sqrt(0.5)*(randn([N_Rx,N_Tx,N_TotUser])+1j*randn([N_Rx,N_Tx,N_TotUser]));


        for iAlg=1:N_Algorithm
            switch iAlg
                case 1
%                     disp('Proposed');
%                     tic
                    SumRate(iNuser,iAlg)=SumRate(iNuser,iAlg)+Proposed(H,N_MaxSupUser,SNR);
%                     toc
                case 2
%                     disp('c-based');
%                     tic
                    SumRate(iNuser,iAlg)=SumRate(iNuser,iAlg)+c_based(H,N_MaxSupUser,SNR);
%                     toc
                case 3
%                     disp('ChodalDistance');
%                     tic
                    SumRate(iNuser,iAlg)=SumRate(iNuser,iAlg)+cd_based(H,N_MaxSupUser,SNR);
%                     toc
                case 4
%                     disp('SRN');
%                     tic
                    SumRate(iNuser,iAlg)=SumRate(iNuser,iAlg)+SRN_based(H,N_MaxSupUser,SNR);
%                     toc
            end
        end

    end
end

SumRate=SumRate/N_Iteration
% SumRate=real(SumRate);

plot(N_User,SumRate(:,1),'-o',N_User,SumRate(:,2),'-*',N_User,SumRate(:,4),'-<',N_User,SumRate(:,3),'->');
legend('Proposed','c-algorithm','psrn-algorithm','cd-algorithm',2);
xlabel('Number of Total User'); ylabel('Average Sum Rate [bps/Hz]'); 
grid on;