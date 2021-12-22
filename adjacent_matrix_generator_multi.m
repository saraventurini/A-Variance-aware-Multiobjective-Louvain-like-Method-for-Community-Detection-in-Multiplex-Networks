%Informative case
%each layer same communities
%SBM with probabilities p and q
%Function to create the adjacent matrix of a undirected graph with m
%communities of size L and k layers
%Input: L array where L(c) is size of community c 
%       (sum(L) matrix size)
%       p probability arc intercommunity
%       q probability arc intracommunity 
%       k number of layers
%Output: M tensor k dim which elements are adjacent matrices of the layers

function [M]=adjacent_matrix_generator_multi(L,p,q,k)
for i=1:k
    M(:,:,i) = adjacent_matrix_generator(L,p,q);
end
end


