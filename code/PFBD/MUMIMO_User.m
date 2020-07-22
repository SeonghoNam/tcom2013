clc; clear all; close all;

addpath('.\Gram-Schmidt Process');
N_Tx=4;
N_Rx=2;

N_Algorithm=5;

N_TotUser=10;
N_MaxSupUser=ceil(N_Tx/N_Rx);


SNR_dB=10;
SNR_real=10.^(SNR_dB*0.1);

N_Iteration=1000;



% N_User=[N_MaxSupUser 4 6 8 10 15 20]; 
N_User=[50 100]; 
SumRate=zeros(length(N_User),N_Algorithm);

for iNuser=1:length(N_User);
    N_TotUser=N_User(iNuser)
    
    TotUserSet=1:N_TotUser;
    
    SNR=SNR_real;
    
    clear TotalCombination;
    for iSupUser=1:N_MaxSupUser
        TotalCombination{iSupUser}=pick(TotUserSet,iSupUser,'');
    end

    
    for iIter=1:N_Iteration

        H=sqrt(0.5)*(randn([N_Rx,N_Tx,N_TotUser])+1j*randn([N_Rx,N_Tx,N_TotUser]));


        for iAlg=1:N_Algorithm
            switch iAlg
                case 1
                    disp('Optimal');
                    tic
                    SumRate(iNuser,iAlg)=SumRate(iNuser,iAlg)+OptimalBD(H,TotalCombination,SNR);
                    toc
                case 2
                    disp('c-based');
                    tic
                    SumRate(iNuser,iAlg)=SumRate(iNuser,iAlg)+c_based(H,N_MaxSupUser,SNR);
                    toc
                case 3
                    disp('ER-based');
                    tic
                    SumRate(iNuser,iAlg)=SumRate(iNuser,iAlg)+ER_based(H,N_MaxSupUser,SNR);
                    toc
                case 4
                    disp('ChodalDistance');
                    tic
                    SumRate(iNuser,iAlg)=SumRate(iNuser,iAlg)+CD_based(H,N_MaxSupUser,SNR);
                    toc
                case 5
                    disp('Angle-based (det)');
                    tic
                    SumRate(iNuser,iAlg)=SumRate(iNuser,iAlg)+CD_based4(H,N_MaxSupUser,SNR);
                    toc
            end
        end

    end
end

SumRate=SumRate/N_Iteration
% SumRate=real(SumRate);

plot(N_User,SumRate,'-o');
legend('Optimal','c-based','ER-based','Distnace-based (CD)','Angle-based (det)');
xlabel('Number of Total User'); ylabel('Sum Capacity [bps/Hz]'); 
