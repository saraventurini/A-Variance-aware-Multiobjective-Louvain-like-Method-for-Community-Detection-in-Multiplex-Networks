%Normalized Mutual Information (NMI)
%INPUT: T confusion matrix
%OUTPUT: n Normalized Mutual Information

function[n]=NMI(T)
%H(C) entropy of C cluster labels
HC=0;
for i=1:size(T,1)
    HC = HC + sum(T(i,:)) * log2(sum(T(i,:))/sum(sum(T)));
end
%H(Y) entropy of class labels
HY=0;
for j=1:size(T,2)
    HY = HY + sum(T(:,j)) * log2(sum(T(:,j))/sum(sum(T)));
end
%I(Y;C) Mutual Information between Y and C
IYC=0;
for i=1:size(T,1)
    for j=1:size(T,2)
        if T(i,j)~= 0
            IYC = IYC + T(i,j) * log2((T(i,j)*sum(sum(T)))/(sum(T(i,:))*sum(T(:,j))));
        end
    end
end
%NMI
n = -(2*IYC)/(HC+HY);
end 