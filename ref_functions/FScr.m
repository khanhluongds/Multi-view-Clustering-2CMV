function [Fscore,RandIdx,Prec,Rec]=FScr(dResults)
%dResults=[1 9;1 9;1 9;1 9;1 9;1 8;2 9;2 8;2 8;2 8;2 8;2 7;3 9;3 9;3 7;3 7;3 7]; 
%Example taken from http://nlp.stanford.edu/IR-book/html/htmledition/evaluation-of-clustering-1.html
%Result should be Precision=0.5, Recall=0.4545, FScore=0.4762, & RandIdx=0.6765
N=size(dResults,1);if N<2, return; end;% dResults=[label, cluster] .. label is the ground truth
beta=1;k=(1+beta^2);% Meaning it calculate F1-Score
TP=0;FN=0;FP=0;%reverseStr='';%initialize all to zeros
for i=1:N-1 %looping through all real class [Ground Truth of the cluster]
    for j=i+1:N %fprintf('i=%d, j=%d \n',i,j);
		%msg=sprintf('%d/%d %d/%d\n',i,N,j,N);fprintf([reverseStr,msg]);
		same_class=dResults(i,1)==dResults(j,1);same_clust=dResults(i,2)==dResults(j,2);
        if same_class&&same_clust
            TP=TP+1;%class and cluster are the same
        elseif same_class&&~same_clust
            FP=FP+1;% assigning two similar documents to different clusters
        elseif ~same_class&&same_clust
            FN=FN+1;%assigning 2 dissimilar documents to the same cluster
        end
        %reverseStr=repmat(sprintf('\b'),1,length(msg));
    end
end%[TP FP TN FN]
TN=(N*(N-1)/2)-(TP+FP+FN);% All Possible Pairwise-(~TN): TN=assigning 2 dissimilar documents to different clusters
if TP+FP==0,Prec=0;else Prec=TP/(TP+FP);end
if TP+FN==0,Rec=0; else Rec=TP/(TP+FN);end
D=(k*TP+(beta^2)*FN+FP);
if D~=0, Fscore=k*TP/D;else Fscore=1;end
RandIdx=(TP+TN)/(TP+FP+FN+TN);%percentage of correct decision: 1 perfect, 0 worse 
