function [pdf x]=plotpdf(data)
    data_linear=reshape(data,[1 numel(data)]);
    [pdf,x]= hist(data_linear,500);  % hist �Լ��� �̿��Ͽ� �������� �׾��ش�

    resol=x(2)-x(1); % resolution�� ����Ѵ�.
    pdf= pdf/(length(data_linear)*resol); % Normalization �� �ش�.
%     cdf=cumsum(pdf);


end
