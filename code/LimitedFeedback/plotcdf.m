function [cdf x]=plotcdf(data)
    data_linear=reshape(data,[1 numel(data)]);
    [pdf,x]= hist(data_linear,1000);  % hist �Լ��� �̿��Ͽ� �������� �׾��ش�

    resol=x(2)-x(1); % resolution�� ����Ѵ�.
    pdf= pdf/(length(data_linear)*resol); % Normalization �� �ش�.
    cdf=cumsum(pdf*resol);


end
