function [idx, quant_val] = myquantiz(val, partitions, codebook)
    unquant_val = val(:);
    my_partitions = reshape(partitions, [1, numel(partitions)]);
    my_codebook = reshape(codebook, [1, numel(codebook)]);
    idx = sum(unquant_val > my_partitions, 2);
    quant_val = my_codebook(idx + 1);
    idx = reshape(idx, size(val));

end
