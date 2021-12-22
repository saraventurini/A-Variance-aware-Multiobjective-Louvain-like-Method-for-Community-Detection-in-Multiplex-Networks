%Function calculates the confusion matrix. It is used in "wrong" function
%INPUT: E vector expected communities
%       P vector predicted communities
%OUTPUT: T confusion matrix 
function [T] = confusion_matrix(E,P)

[E]=reindex_com(E); %reindex community in E 
[P]=reindex_com(P); %reindex community in P 

%columns expected communities - rows predicted communities
Es=sort(E);
Eu=Es([true;diff(Es(:))>0]);
Ps=sort(P);
Pu=Ps([true;diff(Ps(:))>0]);
T = zeros(length(Pu),length(Eu)); 

for i = Eu
    Pi = P(E==i);
    for j = Pu
        T(j,i) = sum(Pi==j);
    end
end

end