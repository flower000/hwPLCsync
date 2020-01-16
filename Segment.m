function [corr] = Segment(seq1,seq2,num)
%Segment: calculate the correlation between the sequence 1 and 2; 
%   
    num = floor(num);
    % zero padding
    if length(seq1) ~= length(seq2)
        corr = 0;
        return;
    else if num ~= 1
        seq1 = [seq1;zeros(num-mod(length(seq1),num),1)];
        seq2 = [seq2;zeros(num-mod(length(seq1),num),1)];
    end
    % reshape
    len = length(seq1);
    temp = reshape(seq1,[num,len/num]) .* conj(reshape(seq2,[num,len/num]));
    % calculate correlation
    if num==1
        corr = sum(abs(temp));
    else
        corr = sum(abs(sum(temp)));
    end
end

