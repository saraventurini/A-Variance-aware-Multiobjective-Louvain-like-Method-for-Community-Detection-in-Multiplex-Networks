%SBM with probabilities p and q
%Function to create the adjacent matrix of a undirected graph with m
%communities of size L 
%Input: L array where L(c) is size of community c 
%       (m*sum(L) matrix size)
%       p probability arc intercommunity
%       q probability arc intracommunity 
%Output: M adjacent matrix of the graph

function [M] = adjacent_matrix_generator(L,p,q)

%M = rand(n) returns an n-by-n matrix of random entries between 0 and 1 
n = sum(L); %dimension matrix
m = length(L); %number communities 
M = rand(n);

%Creation of the communities putting an edge between the nodes with a
%probability p
sum_m=0;
for h=1:m
    if h~=1
       sum_m = sum_m + L(h-1); 
    end
    for i=1:L(h)
        for j=(sum_m+i):(sum_m+L(h))
            if M(i+sum_m,j)<p
                M(i+sum_m,j)=1;
            else
                M(i+sum_m,j)=0;
            end
        end
    end
end

%Creation of the edges between the nodes of different communities
%with a probability q
sum_m=0;
for h=1:(m-1)
    if h~=1
       sum_m = sum_m + L(h-1); 
    end
    for i= (sum_m+1) : (sum_m+L(h)) 
        for j=(sum_m+L(h)+1):n
            if M(i,j)<q
                M(i,j)=1;
            else
                M(i,j)=0;
            end
        end
    %end
    end
end

M = M - diag(diag(M));
%Symetrize matrix 
M=triu(M)+(triu(M,1))';

end
