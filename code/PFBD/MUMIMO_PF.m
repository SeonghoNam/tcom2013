clc; clear all; close all;

addpath('.\Gram-Schmidt Process');
N_Tx=8;
N_Rx=2;

N_Algorithm=5;

N_TotUser=40;
N_MaxSupUser=ceil(N_Tx/N_Rx);

TotUserSet=1:N_TotUser;

PF_Window=100;
PF_Iter=0;

SNR_dB = linspace(0,20,N_TotUser);
SNR_real=10.^(SNR_dB*0.1);
% SNR_real=SNR_real/(sum(SNR_real))*N_TotUser;


N_Iteration=10000;

for iSupUser=1:N_MaxSupUser
    TotalCombination{iSupUser}=pick(TotUserSet,iSupUser,'');
end

RateSet=zeros(N_TotUser,N_Algorithm);
InstaneousRateSet=zeros(N_TotUser,N_Algorithm);
AvgRateSet=ones(N_TotUser,N_Algorithm);



for iIter=1:N_Iteration
    iIter

    H=sqrt(0.5)*(randn([N_Rx,N_Tx,N_TotUser])+1j*randn([N_Rx,N_Tx,N_TotUser]));
    
    for iUser=TotUserSet
        H(:,:,iUser)=sqrt(SNR_real(iUser))*H(:,:,iUser);
    end
    
    for iAlg=1:N_Algorithm
        switch iAlg
            case 1
%                 disp('Sum Rate');
%                 tic
                AvgRateSet(:,iAlg)=1;
                InstaneousRateSet(:,iAlg)=c_based(H,N_MaxSupUser,AvgRateSet(:,iAlg));
%                 toc
            case 2
%                     disp('c_based');
%                     tic
                InstaneousRateSet(:,iAlg)=c_based(H,N_MaxSupUser,AvgRateSet(:,iAlg));
%                     toc
            case 3
%                     disp('prop');
%                     tic
                InstaneousRateSet(:,iAlg)=Proposed_SimpleDet(H,N_MaxSupUser,AvgRateSet(:,iAlg));
%                     toc
            case 4
%                     disp('prop');
%                     tic
                InstaneousRateSet(:,iAlg)=SRN_based(H,N_MaxSupUser,AvgRateSet(:,iAlg));
%                     toc
            case 5
%                     disp('prop');
%                     tic
                InstaneousRateSet(:,iAlg)=Proposed_SimpleDet2(H,N_MaxSupUser,AvgRateSet(:,iAlg));
%                     toc
        end
    end
    
%     if iIter==1
%         AvgRateSet=AvgRateSet+InstaneousRateSet;
%     elseif iIter<PF_Window
%         AvgRateSet=(1-1/iIter)*AvgRateSet+(1-1/iIter)*InstaneousRateSet;
%     else
        AvgRateSet=(1-1/PF_Window)*AvgRateSet+(1-1/PF_Window)*InstaneousRateSet;
%     end

    if iIter>PF_Iter
        RateSet=RateSet+InstaneousRateSet;
    end

end


RateSet=RateSet/(N_Iteration-PF_Iter)

SumRate=sum(RateSet)

FairnessIndex=sum(RateSet).^2./(N_TotUser*sum(RateSet.^2))


plot(TotUserSet,RateSet(:,3),'-o',TotUserSet,RateSet(:,2),'-o',TotUserSet,RateSet(:,4),'-o',TotUserSet,RateSet(:,1),'-ko',TotUserSet,RateSet(:,5),'-d');
legend('Proposed PF','c-algorithm PF','psrn-algorithm PF','SumRate Max','individual c-algorithm',2);
xlabel('User Index'); ylabel('Average Rate [bps/Hz]'); 
ylim([0 1])
