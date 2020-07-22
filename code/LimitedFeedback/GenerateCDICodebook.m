function [Codebook]=GenerateCDICodebook(nt,nr,Nbit)

Qbit=Nbit;
Codebook=cell(1,2^Qbit);
for i=1:2^Qbit
    H=sqrt(1/2)*(randn(nr,nt)+1j*randn(nr,nt));
    [u s v]=svd(H);
    Codebook{i}=v(:,1:nr)';
end



