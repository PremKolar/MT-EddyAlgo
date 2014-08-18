function sumOfAandB=multiDnansum(A,B)
    apnd1d = @(x) reshape(x,[1,size(x)])  ;
    sumOfAandB=squeeze(nanmean([apnd1d(A); apnd1d(B)]));
end