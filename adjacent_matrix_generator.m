%SBM with probabilities p and q
%Function to create the adjacent matrix of a undirected graph with m
%communities of size L 
%Input: m number of community
%       L dimension of each community (every community same size)
%       (m*L matrix size)
%       p probability arc intercommunity
%       q probability arc intracommunity 
%Output: M adjacent matrix of the graph

function [M] = adjacent_matrix_generator(m,L,p,q)

%M = rand(n) returns an n-by-n matrix of random entries between 0 and 1 
n = m*L; %dimension matrix
M = rand(n);

%Creation of the communities putting an edge between the nodes with a
%probability p
for h=0:(m-1)
    for i=1:L
        for j=i:L*(h+1)
            if M(i+L*h,j)<p
                M(i+L*h,j)=1;
            else
                M(i+L*h,j)=0;
            end
        end
    end
end

%Creation of the edges between the nodes of different communities
%with a probability q
for h=0:(m-1)
    for i=1:(L*h) 
        for j=(1+h*L):((h+1)*L)
            if M(i,j)<q
                M(i,j)=1;
            else
                M(i,j)=0;
            end
        end
    %end
    end
end

%Symetrize matrix 
M=triu(M)+(triu(M,1))';

end
