%Informative case
%each layer same communities
%SBM with probabilities p and q
%Function to create the adjacent matrix of a undirected graph with m
%communities of size L and k layers
%Input: m number of community
%       L dimension of each community (every community same size)
%       (m*L each layer matrix size)
%       p probability arc intercommunity
%       q probability arc intracommunity 
%       k number of layers
%Output: M tensor k dim which elements are adjacent matrices of the layers

function [M]=adjacent_matrix_generator_multi(m,L,p,q,k)
for i=1:k
    M(:,:,i) = adjacent_matrix_generator(m,L,p,q);
end
end


