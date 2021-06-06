%Noisy case
%First layer random (p==q) and all the other layers informative p>>q 
%Input: m number of community
%       L dimension of each community (every community same size)
%       (m*L each layer matrix size)
%       p probability arc intercommunity (informative layers)
%       q probability arc intracommunity (informative layers)
%       s probability arc inter and intra community (noise layers)
%       k number of total layers
%Output: M tensor k dim which elements are adjacent matrices of the layers

function [M]=adjacent_matrix_generator_multi_r(m,L,p,q,s,k)

M(:,:,1) = adjacent_matrix_generator(m,L,s,s); %noise

if k>1
for i=2:k
    M(:,:,i) = adjacent_matrix_generator(m,L,p,q); %informative
end
end
end