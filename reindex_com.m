%function to reindex communities. It is used in "confusion_matrix" function
function [C] = reindex_com(COMold) 

C = sparse(1,length(COMold)); %vector with a index for each node 

%unique communities
COMolds=sort(COMold);
COMu=COMolds([true;diff(COMolds(:))>0]);

for l=1:length(COMu)
    C(COMold==COMu(l)) = l; %Riname the communities  
end
end