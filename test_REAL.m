% Tests on Real Datasets
function test_REAL
oldFolder = cd;
dataset=["3sources.mat", "BBCSport2view_544.mat","cora.mat","UCI_mfeat.mat","WikipediaArticles.mat","aucs.mat","dkpol.mat"];

for data = dataset
cd ./Real_datasets
load(data); %load data 
cd(oldFolder)
if data == "aucs.mat" %node_nolabels
    %informative case
    REAL(W_cell,labels,data,node_nolabels) 
    %noisy case 2 for s = 0.01, 0.03, 0.05
    REAL_n2(W_cell,labels,data,0.01,node_nolabels); 
    REAL_n2(W_cell,labels,data,0.03,node_nolabels); 
    REAL_n2(W_cell,labels,data,0.05,node_nolabels); 
    %noisy case 3 for s = 0.01, 0.03, 0.05
    REAL_n3(W_cell,labels,data,0.01,node_nolabels); 
    REAL_n3(W_cell,labels,data,0.03,node_nolabels); 
    REAL_n3(W_cell,labels,data,0.05,node_nolabels); 
else
    %informative case
    REAL(W_cell,labels,data) 
    %noisy case 2 for s = 0.01, 0.03, 0.05
    REAL_n2(W_cell,labels,data,0.01); 
    REAL_n2(W_cell,labels,data,0.03); 
    REAL_n2(W_cell,labels,data,0.05); 
    %noisy case 3 for s = 0.01, 0.03, 0.05
    REAL_n3(W_cell,labels,data,0.01); 
    REAL_n3(W_cell,labels,data,0.03); 
    REAL_n3(W_cell,labels,data,0.05); 
end
end
cd(oldFolder)
end
