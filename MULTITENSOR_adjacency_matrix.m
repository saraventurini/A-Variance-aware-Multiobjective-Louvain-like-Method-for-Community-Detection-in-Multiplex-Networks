% From adjecency matrix to file for MULTITENSOR method 
% `E node1 node2 3 0 0 1`
% The first columns tells the algorithm that the row denotes an edge; the second and third are the source and target nodes of that edge, respectively; l+3 column tells if there is that edge in the l-th layer and the weigth (must be integer). 
% In this example the edge node1 --> node2 exists in layer 1 with weight 3 and in layer 4 with weight 1, but not in layer 2 and 3.
% INPUT: M adjacency matrix tensor. A(:,:,l) is adjacency matrix l-layer.
%        n number of nodes. 


function MULTITENSOR_adjacency_matrix(M,n)
oldFolder = cd;
cd ./MultiTensor-master/data 
%{
%REALDATASETS
%download real dataset
%W_cell cell s.t. W_cell{k} adjacency matrix layer k
%labels array 1xN true communities s.t. labels(n) community node n 
%load('C:\Users\Sara\Desktop\Codes\Real_datasets\3sources.mat','labels','W_cell');
load('C:\Users\Sara\Desktop\Codes\Real_datasets\3sources.mat','labels','W_cell');
% datasets = ['3sources.mat', 'BBCSport2view_544.mat', 'cora.mat', 'UCI_mfeat.mat', 'WikipediaArticles.mat'];
%W_cell, labels

n = size(W_cell{1},1); %numbe rof nodes
K = size(W_cell,1); %number of layers
for k=1:K
    M(:,:,k)=full(W_cell{k}); %Matrix
    M(:,:,k)=M(:,:,k)-diag(diag(M(:,:,k))); % no loop
end
M_sum = sum(M(:,:,:),3); %summ M tensor 3-dim 
%}

%{
%EXAMPLE
n = 2;
K = 2;
for k=1:K
    M(:,:,k)=ones(n); %Matrix
end
M(:,:,3)=eye(n);
M_sum = sum(M(:,:,:),3); %summ M tensor 3-dim 
%}

M_sum = sum(M(:,:,:),3); %summ M tensor 3-dim
% Write edjes on file
e = 0;%number edges
for i=1:n
    for j=1:n
        if M_sum(i,j)~=0
           e = e+1;
           Edjes_matrix(e,:) = cat(2,[i,j],reshape(M(i,j,:), 1, [])); 
        end
    end
end

writematrix(Edjes_matrix,'adjacencyD.txt','Delimiter',' ');

%'E' at beginning each line
fid1 = fopen('adjacencyD.txt');
fid2 = fopen('adjacency.txt','w');
tline = fgets(fid1);
Character2add = 'E'; % here put the character to be added
while ischar(tline)
    %disp(tline)
    fprintf(fid2,'%s %s',Character2add,tline);
    tline = fgets(fid1);
end
fclose(fid1);
fclose(fid2);
delete('adjacencyD.txt')

%.dat file 
[~, baseFileName, ~] = fileparts('adjacency.txt');
newFullFileName = sprintf('%s.dat', baseFileName);
copyfile('adjacency.txt',newFullFileName)
delete('adjacency.txt')
cd(oldFolder)
end
