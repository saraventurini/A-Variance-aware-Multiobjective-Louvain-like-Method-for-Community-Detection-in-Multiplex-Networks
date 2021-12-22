%SBM Noisy case
%First layer random (p==q) and all the other layers informative p>>q 
%Input: L array where L(c) is size of community c
%       (sum(L) each layer matrix size)
%       p probability arc intercommunity (informative layers)
%       q probability arc intracommunity (informative layers)
%       s probability arc inter and intra community (noise layers)
%       k1 number of informative layers
%       k2 number of noisy layers 
%       (k1+k2 total number of layers)
%Output: M tensor k dim which elements are adjacent matrices of the layers

function [M]=adjacent_matrix_generator_multi_r(L,p,q,s,k1,k2)

%informative
for i=1:k1
    M(:,:,i) = adjacent_matrix_generator(L,p,q); 
end
%noise 
for i=(k1+1):(k1+k2)
    M(:,:,i) = adjacent_matrix_generator(L,s,s); 
end
end



