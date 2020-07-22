function [QuantizedValue] = CDIQuantizer(sample,Codebook)


nr=size(sample,1);
[u s v]=svd(sample);
H=v(:,1:nr)';
for i=1:length(Codebook)
    cd(i)=sqrt(nr-trace(H*Codebook{i}'*Codebook{i}*H'));
end
[temp index]=min(cd);
QuantizedValue=Codebook{index};
end

