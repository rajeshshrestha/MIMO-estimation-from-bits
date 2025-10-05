function [codes] = generateOVSF(Nt, N_pilots)
%GENERATEOVSF Generate 'N_pilots' OVSF codes each of length 'code_size'

% N_pilots should be less than or equal to 2^Nt
assert(N_pilots >= Nt, 'Invalid combination of Nt and code N_pilots');

% check if N_pilots is a power of 2
assert(log2(N_pilots) == round(log2(N_pilots)), 'N_pilots should be a power of 2');

Y = [[1]];
for level=1:log2(N_pilots)
    X=Y;
    Y=[];
    for i=1:size(X,1)
        Y=[Y; [X(i,:) X(i,:)]; [X(i,:) -X(i,:)]];
    end
end

s = RandStream('mlfg6331_64');
codes = datasample(s, Y, Nt, 1, 'Replace', false);
codes=codes(:,1:N_pilots)/sqrt(N_pilots); % Normalize the codes 
end
