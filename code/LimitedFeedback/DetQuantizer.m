function [QuantizedValue] = DetQuantizer(sample,level,value)

for i=1:length(value)
    if sample > level(i) && sample < level(i+1)
        QuantizedValue=value(i);
        continue;
    end
end

end

