function out=calc_information_entropy(I)
[N,~] = histcounts(I(:),255);
N=N/sum(N(:));
p=-1*N.*log(N);
p(isnan(p))=0;
out=sum(p(:));
end