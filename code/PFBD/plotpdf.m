function [pdf x]=plotpdf(data)
    data_linear=reshape(data,[1 numel(data)]);
    [pdf,x]= hist(data_linear,500);  % hist 함수를 이용하여 구간별로 쌓아준다

    resol=x(2)-x(1); % resolution을 계산한다.
    pdf= pdf/(length(data_linear)*resol); % Normalization 해 준다.
%     cdf=cumsum(pdf);


end
