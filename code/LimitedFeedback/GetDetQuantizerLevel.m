function [level value]=GetDetQuantizerLevel(SNR,nt,nr,Nbit)

Qbit=Nbit;
for i=1:10000
    H=sqrt(SNR/2)*(randn(nr,nt)+1j*randn(nr,nt));
    data(i)=real(det(H*H'));
end

[y x]=plotcdf(data);

b(1)=0;
b(2^Qbit+1)=1;

for i=1:2^Qbit
    b(i+1)=i/2^Qbit;
    [temp index]=min(abs(y-(b(i+1))));
    a(i+1)=x(index);
    value(i)=(a(i+1)+a(i))/2;

end
a(1)=0;
a(2^Qbit+1)=inf;
% value(2^Qbit)=(40+a(2^Qbit))/2;
if Nbit==0
    level=[0 inf];
    value=1;
end
level=a;

