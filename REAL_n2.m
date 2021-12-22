%Real Datasets-Noisy Case: all informative layers + one noisy layer
%The function tests each methods 10 times becouse the noise is added in a
%stochastic way. 
%Each method is tested 10 times on the same graph (becouse in the noisy case the communities are reindexed randomly). 
%It prints the average results in a csv file called "data.csv"
%INPUT: W_cell cell containing adjacent matrix of each layer
%       labels vector right communities 
%       data name datatset
%       s value of p_in=p_out for the random layer
%       node_nolabels (just for aucs dataset) list of nodes with no true community assigned, then no consider for accuracy calculation
%Choose between lines 23 and 24 depending on the dataset

function REAL_n2(W_cell,labels,data,s,node_nolabels)

Methods={'Time';'Acc';'NMI'};

%Martrix
N = size(W_cell{1},1); %nodes
if data=="aucs.mat" %if nargin==4
    n = N-length(node_nolabels); %no consider nodes without label in calculation of accuracy
else
    n = N;
end
k = size(W_cell,1)+1; %layers
%informative layers
for i=1:(k-1)
    M(:,:,i)=full(W_cell{i});
    A{i}=M(:,:,i);
end


%Expected comunities
if data=="aucs.mat" %if nargin==4
    labels(node_nolabels)=[]; %delate community 0 of nodes without label
end
c = labels; 

if any(strcmp(["3sources.mat"; "BBCSport2view_544.mat" ; "WikipediaArticles.mat"; "aucs.mat" ; "dkpol.mat"],data))
m = size(unique(labels),2); %use for 3sources, BBC,BBCSport, Wikipedia
else
m = size(unique(labels),1); %use for UCI, cora
end

%state-of-the-art methods
%CoReg
addpath ./utils
addpath ./MinMaxSelection
addpath ./CoReg
lambda=1e-2;
lambda = ones(k, 1) * lambda;
%AWP
addpath ./AWP
opts.NITER = 300;
%MCGC
addpath ./MCGC
beta=0.6;
%PM
addpath(genpath('utils'))
addpath(genpath('subroutines'))
addpath(genpath('ThePowerMeanLaplacianForMultilayerGraphClustering'))
z = -10;
%MT
addpath ./MultiTensor-master/python
addpath ./MultiTensor-master/data
%SCML
addpath ./scml
lambda2 = 0.5;
%PMM
addpath ./PMM
ell = 10; %ell = fix(n/2);
maxKmeans = 5;

Time1_m=[];
Acc1_m=[];
nmi1_m=[];
Time2_m=[];
Acc2_m=[];
nmi2_m=[];
Time3_m=[];
Acc3_m=[];
nmi3_m=[];
Time4_m=[];
Acc4_m=[];
nmi4_m=[];
Time5_m=[];
Acc5_m=[];
nmi5_m=[];
Time6_m=[];
Acc6_m=[];
nmi6_m=[];
Time7_m=[];
Acc7_m=[];
nmi7_m=[];
Time8_m=[];
Acc8_m=[];
nmi8_m=[];
Time9_m=[];
Acc9_m=[];
nmi9_m=[];
Time10_m=[];
Acc10_m=[];
nmi10_m=[];
Time11_m=[];
Acc11_m=[];
nmi11_m=[];
Time12_m=[];
Acc12_m=[];
nmi12_m=[];
Time13_m=[];
Acc13_m=[];
nmi13_m=[];
Time14_m=[];
Acc14_m=[];
nmi14_m=[];
Time15_m=[];
Acc15_m=[];
nmi15_m=[];
Time16_m=[];
Acc16_m=[];
nmi16_m=[];
Time17_m=[];
Acc17_m=[];
nmi17_m=[];
Time18_m=[];
Acc18_m=[];
nmi18_m=[];
Time19_m=[];
Acc19_m=[];
nmi19_m=[];
Time20_m=[];
Acc20_m=[];
nmi20_m=[];
Time21_m=[];
Acc21_m=[];
nmi21_m=[];
Time22_m=[];
Acc22_m=[];
nmi22_m=[];
Time23_m=[];
Acc23_m=[];
nmi23_m=[];
Time24_m=[];
Acc24_m=[];
nmi24_m=[];
Time25_m=[];
Acc25_m=[];
nmi25_m=[];
Time26_m=[];
Acc26_m=[];
nmi26_m=[];

for N_iter2 = 1:10

Time1=[];
Acc1=[];
nmi1=[];
Time2=[];
Acc2=[];
nmi2=[];
Time3=[];
Acc3=[];
nmi3=[];
Time4=[];
Acc4=[];
nmi4=[];
Time5=[];
Acc5=[];
nmi5=[];
Time6=[];
Acc6=[];
nmi6=[];
Time7=[];
Acc7=[];
nmi7=[];
Time8=[];
Acc8=[];
nmi8=[];
Time9=[];
Acc9=[];
nmi9=[];
Time10=[];
Acc10=[];
nmi10=[];
Time11=[];
Acc11=[];
nmi11=[];
Time12=[];
Acc12=[];
nmi12=[];
Time13=[];
Acc13=[];
nmi13=[];
Time14=[];
Acc14=[];
nmi14=[];
Time15=[];
Acc15=[];
nmi15=[];
Time16=[];
Acc16=[];
nmi16=[];
Time17=[];
Acc17=[];
nmi17=[];
Time18=[];
Acc18=[];
nmi18=[];
Time19=[];
Acc19=[];
nmi19=[];
Time20=[];
Acc20=[];
nmi20=[];
Time21=[];
Acc21=[];
nmi21=[];
Time22=[];
Acc22=[];
nmi22=[];
Time23=[];
Acc23=[];
nmi23=[];
Time24=[];
Acc24=[];
nmi24=[];
Time25=[];
Acc25=[];
nmi25=[];
Time26=[];
Acc26=[];
nmi26=[];

%add one noise layer
M(:,:,k) = adjacent_matrix_generator([N],s,s);
A{k}=M(:,:,k);
%file adjacency matrix for MultiTensor method
MULTITENSOR_adjacency_matrix(M,N);
%file adjacency matrix for INFOMAP method "adjacency_INFOMAP.net"
INFOMAP_adjacency_matrix(M,N,k);

%Run the methods on the matrix
%Generalized Louvain
for N_iter = 1:5 
            tx1 = tic;
    [COMTY1, ending1] = GL_r(M); 
    tEnd1 = toc(tx1);
    if ending1==0
    Time1 = [Time1 tEnd1]; %time
    if data=="aucs.mat" %if nargin==4
        COMTY1.COM{end}(node_nolabels)=[]; %no consider nodes without label in calculation of accuracy
    end 
    Acc1 = [Acc1 (n-wrong(c,COMTY1.COM{end}))/(n)]; %accuracy 
    [T1] = confusion_matrix(c,COMTY1.COM{end}'); %confusion matrix
    nmi1 = [nmi1 NMI(T1)]; %NMI
    end
end
Time1_m = [Time1_m mean(Time1)];
Acc1_m = [Acc1_m mean(Acc1) ];
nmi1_m = [nmi1_m mean(nmi1)];

%Louvain Expansion Method Functin F+, lambda=0.1
for N_iter = 1:5
            tx2 = tic;
    [COMTY2, ending2] = EVP(M,0.1);
    tEnd2 = toc(tx2);
    if ending2==0
    Time2 = [Time2 tEnd2];
    if data=="aucs.mat" %if nargin==4
        COMTY2.COM{end}(node_nolabels)=[]; %no consider nodes without label in calculation of accuracy
    end 
    Acc2 = [Acc2 (n-wrong(c,COMTY2.COM{end}))/(n)];
    [T2] = confusion_matrix(c,COMTY2.COM{end}');
    nmi2 = [nmi2 NMI(T2)];
    end
end
Time2_m = [Time2_m mean(Time2)];
Acc2_m = [Acc2_m mean(Acc2) ];
nmi2_m = [nmi2_m mean(nmi2)];
%Louvain Expansion Method Functin F+, lambda=0.3
for N_iter = 1:5
        tx3 = tic;
    [COMTY3, ending3] = EVP(M,0.3); 
    tEnd3 = toc(tx3);
    if ending3==0
    Time3 = [Time3 tEnd3];
    if data=="aucs.mat" %if nargin==4
        COMTY3.COM{end}(node_nolabels)=[]; %no consider nodes without label in calculation of accuracy
    end 
    Acc3 = [Acc3 (n-wrong(c,COMTY3.COM{end}))/(n)];
    [T3] = confusion_matrix(c,COMTY3.COM{end}');
    nmi3 = [nmi3 NMI(T3)];
    end
end
Time3_m = [Time3_m mean(Time3)];
Acc3_m = [Acc3_m mean(Acc3) ];
nmi3_m = [nmi3_m mean(nmi3)];
%Louvain Expansion Method Functin F+, lambda=0.5   
for N_iter = 1:5
        tx4 = tic;
    [COMTY4, ending4] = EVP(M,0.5);  
    tEnd4 = toc(tx4);
    if ending4==0
    Time4 = [Time4 tEnd4];
    if data=="aucs.mat" %if nargin==4
        COMTY4.COM{end}(node_nolabels)=[]; %no consider nodes without label in calculation of accuracy
    end 
    Acc4 = [Acc4 (n-wrong(c,COMTY4.COM{end}))/(n)];
    [T4] = confusion_matrix(c,COMTY4.COM{end}');
    nmi4 = [nmi4 NMI(T4)];
    end
end
Time4_m = [Time4_m mean(Time4)];
Acc4_m = [Acc4_m mean(Acc4) ];
nmi4_m = [nmi4_m mean(nmi4)];
%Louvain Expansion Method Functin F+, lambda=0.7
for N_iter = 1:5
        tx5 = tic;
    [COMTY5, ending5] = EVP(M,0.7);  
    tEnd5 = toc(tx5);
    if ending5==0
    Time5 = [Time5 tEnd5];
    if data=="aucs.mat" %if nargin==4
        COMTY5.COM{end}(node_nolabels)=[]; %no consider nodes without label in calculation of accuracy
    end
    Acc5 = [Acc5 (n-wrong(c,COMTY5.COM{end}))/(n)];
    [T5] = confusion_matrix(c,COMTY5.COM{end}');
    nmi5 = [nmi5 NMI(T5)];
    end
end
Time5_m = [Time5_m mean(Time5)];
Acc5_m = [Acc5_m mean(Acc5) ];
nmi5_m = [nmi5_m mean(nmi5)];
%Louvain Expansion Method Functin F+, lambda=0.9
for N_iter = 1:5
        tx6 = tic;
    [COMTY6, ending6] = EVP(M,0.9);  
    tEnd6 = toc(tx6);
    if ending6==0
    Time6 = [Time6 tEnd6];
    if data=="aucs.mat" %if nargin==4
        COMTY6.COM{end}(node_nolabels)=[]; %no consider nodes without label in calculation of accuracy
    end 
    Acc6 = [Acc6 (n-wrong(c,COMTY6.COM{end}))/(n)];
    [T6] = confusion_matrix(c,COMTY6.COM{end}');
    nmi6 = [nmi6 NMI(T6)];
    end
end
Time6_m = [Time6_m mean(Time6)];
Acc6_m = [Acc6_m mean(Acc6) ];
nmi6_m = [nmi6_m mean(nmi6)];
%Louvain Multiobjective Method Average, h=2
for N_iter = 1:5
        tx7 = tic;
    [L7, ending7] = MA_r(M,2); 
    tEnd7 = toc(tx7);
    if ending7
    Time7 = [Time7 tEnd7];
    if data=="aucs.mat" %if nargin==4
        L7{1}{2}{2}(node_nolabels)=[]; %no consider nodes without label in calculation of accuracy
    end 
    Acc7 = [Acc7 (n-wrong(c,L7{1}{2}{2}'))/(n)];
    [T7] = confusion_matrix(c,L7{1}{2}{2});
    nmi7 = [nmi7 NMI(T7)];
    end
end
Time7_m = [Time7_m mean(Time7)];
Acc7_m = [Acc7_m mean(Acc7) ];
nmi7_m = [nmi7_m mean(nmi7)];
%Louvain Multiobjective Method Average, h=3
for N_iter = 1:5
    tx8 = tic;
    [L8, ending8] = MA_r(M,3); 
    tEnd8 = toc(tx8);
    if ending8
    Time8 = [Time8 tEnd8];
    if data=="aucs.mat" %if nargin==4
        L8{1}{2}{2}(node_nolabels)=[]; %no consider nodes without label in calculation of accuracy
    end
    Acc8 = [Acc8 (n-wrong(c,L8{1}{2}{2}'))/(n)];
    [T8] = confusion_matrix(c,L8{1}{2}{2});
    nmi8 = [nmi8 NMI(T8)];
    end
end
Time8_m = [Time8_m mean(Time8)];
Acc8_m = [Acc8_m mean(Acc8) ];
nmi8_m = [nmi8_m mean(nmi8)];    
%Louvain Multiobjective Method Function F+, lambda=0.1, h=2    
for N_iter = 1:5
    tx9 = tic;
    [L9, ending9] = MVP(M,0.1,2); 
    tEnd9 = toc(tx9);
    if ending9
    Time9 = [Time9 tEnd9];
    if data=="aucs.mat" %if nargin==4
        L9{1}{2}{2}(node_nolabels)=[]; %no consider nodes without label in calculation of accuracy
    end
    Acc9 = [Acc9 (n-wrong(c,L9{1}{2}{2}'))/(n)];
    [T9] = confusion_matrix(c,L9{1}{2}{2});
    nmi9 = [nmi9 NMI(T9)];   
    end 
end
Time9_m = [Time9_m mean(Time9)];
Acc9_m = [Acc9_m mean(Acc9) ];
nmi9_m = [nmi9_m mean(nmi9)];
%Louvain Multiobjective Method Function F+, lambda=0.1, h=3    
for N_iter = 1:5
    tx10 = tic;
    [L10, ending10] = MVP(M,0.1,3); 
    tEnd10 = toc(tx10);
    if ending10
    Time10 = [Time10 tEnd10];
    if data=="aucs.mat" %if nargin==4
        L10{1}{2}{2}(node_nolabels)=[]; %no consider nodes without label in calculation of accuracy
    end
    Acc10 = [Acc10 (n-wrong(c,L10{1}{2}{2}'))/(n)];
    [T10] = confusion_matrix(c,L10{1}{2}{2});
    nmi10 = [nmi10 NMI(T10)];   
    end
end
Time10_m = [Time10_m mean(Time10)];
Acc10_m = [Acc10_m mean(Acc10) ];
nmi10_m = [nmi10_m mean(nmi10)];
%Louvain Multiobjective Method Function F+, lambda=0.3, h=2    
for N_iter = 1:5
    tx11 = tic;
    [L11, ending11] = MVP(M,0.3,2); 
    tEnd11 = toc(tx11);
    if ending11
    Time11 = [Time11 tEnd11];
    if data=="aucs.mat" %if nargin==4
        L11{1}{2}{2}(node_nolabels)=[]; %no consider nodes without label in calculation of accuracy
    end
    Acc11 = [Acc11 (n-wrong(c,L11{1}{2}{2}'))/(n)];
    [T11] = confusion_matrix(c,L11{1}{2}{2});
    nmi11 = [nmi11 NMI(T11)];  
    end
end
Time11_m = [Time11_m mean(Time11)];
Acc11_m = [Acc11_m mean(Acc11) ];
nmi11_m = [nmi11_m mean(nmi11)];
%Louvain Multiobjective Method Function F+, lambda=0.3, h=3    
for N_iter = 1:5
    tx12 = tic;
    [L12, ending12] = MVP(M,0.3,3); 
    tEnd12 = toc(tx12);
    if data=="aucs.mat" %if nargin==4
        L12{1}{2}{2}(node_nolabels)=[]; %no consider nodes without label in calculation of accuracy
    end
    if ending12
    Time12 = [Time12 tEnd12];
    Acc12 = [Acc12 (n-wrong(c,L12{1}{2}{2}'))/(n)];
    [T12] = confusion_matrix(c,L12{1}{2}{2});
    nmi12 = [nmi12 NMI(T12)]; 
    end
end
Time12_m = [Time12_m mean(Time12)];
Acc12_m = [Acc12_m mean(Acc12) ];
nmi12_m = [nmi12_m mean(nmi12)];
%Louvain Multiobjective Method Function F+, lambda=0.5, h=2    
for N_iter = 1:5
    tx13 = tic;
    [L13, ending13] = MVP(M,0.5,2); 
    tEnd13 = toc(tx13);
    if ending13
    Time13 = [Time13 tEnd13];
    if data=="aucs.mat" %if nargin==4
        L13{1}{2}{2}(node_nolabels)=[]; %no consider nodes without label in calculation of accuracy
    end
    Acc13 = [Acc13 (n-wrong(c,L13{1}{2}{2}'))/(n)];
    [T13] = confusion_matrix(c,L13{1}{2}{2});
    nmi13 = [nmi13 NMI(T13)]; 
    end
end
Time13_m = [Time13_m mean(Time13)];
Acc13_m = [Acc13_m mean(Acc13) ];
nmi13_m = [nmi13_m mean(nmi13)];
%Louvain Multiobjective Method Function F+, lambda=0.5, h=3    
for N_iter = 1:5
    tx14 = tic;
    [L14, ending14] = MVP(M,0.5,3); 
    tEnd14 = toc(tx14);
    if ending14
    Time14 = [Time14 tEnd14];
    if data=="aucs.mat" %if nargin==4
        L14{1}{2}{2}(node_nolabels)=[]; %no consider nodes without label in calculation of accuracy
    end
    Acc14 = [Acc14 (n-wrong(c,L14{1}{2}{2}'))/(n)];
    [T14] = confusion_matrix(c,L14{1}{2}{2});
    nmi14 = [nmi14 NMI(T14)];   
    end
end
Time14_m = [Time14_m mean(Time14)];
Acc14_m = [Acc14_m mean(Acc14)];
nmi14_m = [nmi14_m mean(nmi14)];
%Louvain Multiobjective Method Function F+, lambda=0.7, h=2
for N_iter = 1:5
    tx15 = tic;
    [L15, ending15] = MVP(M,0.7,2); 
    tEnd15 = toc(tx15);
    if ending15
    Time15 = [Time15 tEnd15];
    if data=="aucs.mat" %if nargin==4
        L15{1}{2}{2}(node_nolabels)=[]; %no consider nodes without label in calculation of accuracy
    end
    Acc15 = [Acc15 (n-wrong(c,L15{1}{2}{2}'))/(n)];
    [T15] = confusion_matrix(c,L15{1}{2}{2});
    nmi15 = [nmi15 NMI(T15)];
    end
end
Time15_m = [Time15_m mean(Time15)];
Acc15_m = [Acc15_m mean(Acc15)];
nmi15_m = [nmi15_m mean(nmi15)];
%Louvain Multiobjective Method Function F+, lambda=0.7, h=3
for N_iter = 1:5
    tx16 = tic;
    [L16, ending16] = MVP(M,0.7,3); 
    tEnd16 = toc(tx16);
    if ending16
    Time16 = [Time16 tEnd16];
    if data=="aucs.mat" %if nargin==4
        L16{1}{2}{2}(node_nolabels)=[]; %no consider nodes without label in calculation of accuracy
    end
    Acc16 = [Acc16 (n-wrong(c,L16{1}{2}{2}'))/(n)];
    [T16] = confusion_matrix(c,L16{1}{2}{2});
    nmi16 = [nmi16 NMI(T16)]; 
    end
end
Time16_m = [Time16_m mean(Time16)];
Acc16_m = [Acc16_m mean(Acc16)];
nmi16_m = [nmi16_m mean(nmi16)];
%Louvain Multiobjective Method Function F+, lambda=0.9, h=2
for N_iter = 1:5
        tx17 = tic;
    [L17, ending17] = MVP(M,0.9,2); 
    tEnd17 = toc(tx17);
    if ending17
    Time17 = [Time17 tEnd17];
    if data=="aucs.mat" %if nargin==4
        L17{1}{2}{2}(node_nolabels)=[]; %no consider nodes without label in calculation of accuracy
    end
    Acc17 = [Acc17 (n-wrong(c,L17{1}{2}{2}'))/(n)];
    [T17] = confusion_matrix(c,L17{1}{2}{2});
    nmi17 = [nmi17 NMI(T17)];
    end
end
Time17_m = [Time17_m mean(Time17)];
Acc17_m = [Acc17_m mean(Acc17)];
nmi17_m = [nmi17_m mean(nmi17)];
%Louvain Multiobjective Method Function F+, lambda=0.9, h=3
for N_iter = 1:5
    tx18 = tic;
    [L18, ending18] = MVP(M,0.9,3); 
    tEnd18 = toc(tx18);
    if data=="aucs.mat" %if nargin==4
        L18{1}{2}{2}(node_nolabels)=[]; %no consider nodes without label in calculation of accuracy
    end
    if ending18
    Time18 = [Time18 tEnd18];
    Acc18 = [Acc18 (n-wrong(c,L18{1}{2}{2}'))/(n)];
    [T18] = confusion_matrix(c,L18{1}{2}{2});
    nmi18 = [nmi18 NMI(T18)];
    end
end
Time18_m = [Time18_m mean(Time18)];
Acc18_m = [Acc18_m mean(Acc18)];
nmi18_m = [nmi18_m mean(nmi18)];
%CoReg
for N_iter = 1:5
    tx19 = tic;
    [label19] = spectral_centroid_multiview_onkernel(A, k, m, lambda);
    tEnd19 = toc(tx19);
    Time19 = [Time19 tEnd19];
    if data=="aucs.mat" %if nargin==4
        label19(node_nolabels)=[]; %no consider nodes without label in calculation of accuracy
    end
    Acc19 = [Acc19 (n-wrong(c,label19'))/(n)];
    [T19] = confusion_matrix(c,label19);
    nmi19 = [nmi19 NMI(T19)];
end
Time19_m = [Time19_m mean(Time19)];
Acc19_m = [Acc19_m mean(Acc19)];
nmi19_m = [nmi19_m mean(nmi19)];
%AWP
for N_iter = 1:5
       tx20 = tic;
    for t=1:k
        embedding20{t} = spectral_embedding(A{t}, m);
    end
    label20 = AWP(embedding20, opts);
    tEnd20 = toc(tx20);
    Time20 = [Time20 tEnd20];
    if data=="aucs.mat" %if nargin==4
        label20(node_nolabels)=[]; %no consider nodes without label in calculation of accuracy
    end
    Acc20 = [Acc20 (n-wrong(c,label20'))/(n)];
    [T20] = confusion_matrix(c,label20);
    nmi20 = [nmi20 NMI(T20)];
end
Time20_m = [Time20_m mean(Time20)];
Acc20_m = [Acc20_m mean(Acc20)];
nmi20_m = [nmi20_m mean(nmi20)];
%MCGC
for N_iter = 1:5
        tx21 = tic;
    [label21] = obj_MVCC(M,k,m,beta);
    tEnd21 = toc(tx21);
    Time21 = [Time21 tEnd21];
    if data=="aucs.mat" %if nargin==4
        label21(node_nolabels)=[]; %no consider nodes without label in calculation of accuracy
    end
    Acc21 = [Acc21 (n-wrong(c,label21'))/(n)];
    [T21] = confusion_matrix(c,label21);
    nmi21 = [nmi21 NMI(T21)]; 
end
Time21_m = [Time21_m mean(Time21)];
Acc21_m = [Acc21_m mean(Acc21)];
nmi21_m = [nmi21_m mean(nmi21)];
%PM
for N_iter = 1:5
        tx22 = tic;
    [label22] = clustering_multilayer_graphs_with_power_mean_laplacian(A, z, m);
    tEnd22 = toc(tx22);
    Time22 = [Time22 tEnd22];
    if data=="aucs.mat" %if nargin==4
        label22(node_nolabels)=[]; %no consider nodes without label in calculation of accuracy
    end
    Acc22 = [Acc22 (n-wrong(c,label22'))/(n)];
    [T22] = confusion_matrix(c,label22');
    nmi22 = [nmi22 NMI(T22)];
end
Time22_m = [Time22_m mean(Time22)];
Acc22_m = [Acc22_m mean(Acc22)];
nmi22_m = [nmi22_m mean(nmi22)];
%MT
%for N_iter = 1:5
    oldFolder = cd;
    cd ./MultiTensor-master/python 
    %run the method    
    python_file_name = ['python main.py -k=' num2str(m) ' -l=' num2str(k) ' -a="adjacency.dat" -E="_endfile.dat"'];    
    tx23 = tic;
    system(python_file_name);
    tEnd23 = toc(tx23);
    %read the result
    cd(oldFolder)
    file_results = ['./MultiTensor-master/data/u_K' num2str(m) '_endfile.dat'];
    B = readtable(file_results);
    B = B{:,:};
    B(:,1)=[];
    [~, label23] = max(B,[],2);
    Time23 = [Time23 tEnd23];
    if data=="aucs.mat" %if nargin==4
        label23(node_nolabels)=[]; %no consider nodes without label in calculation of accuracy
    end
    Acc23 = [Acc23 (n-wrong(c,label23'))/(n)];
    [T23] = confusion_matrix(c,label23');
    nmi23 = [nmi23 NMI(T23)];
%end
Time23_m = [Time23_m mean(Time23)];
Acc23_m = [Acc23_m mean(Acc23)];
nmi23_m = [nmi23_m mean(nmi23)];
%INFOMAP
%for N_iter = 1:5
     oldFolder = cd;
    cd ./INFOMAP
    tx24 = tic;
    python_file_name = 'INFOMAP_method.py';     
    system(python_file_name);
    Time24 = toc(tx24);
    %read the result
    file_results = 'INFOMAP_result.csv';
    B = readtable(file_results);  
    B = B{:,:}; 
    B(:,1)=[]; %node_id layer_id module_id
    layers = B(:,2);
    subB = B(layers==1,[1 3]); %submatrix s layer
    [~,idx] = sort(subB(:,1)); %sort nodes
    sortedmat = subB(idx,:); %sort the submatrix 
    labels_layers = sortedmat(:,2)'; 
    cd(oldFolder)
    [label24] = labels_layers;
    if data=="aucs.mat" %if nargin==4
        label24(node_nolabels)=[]; %no consider nodes without label in calculation of accuracy
    end
    Acc24 = (n-wrong(c,label24))/(n);
    [T24] = confusion_matrix(c,label24);
    nmi24 = NMI(T24);
%end
Time24_m = [Time24_m mean(Time24)];
Acc24_m = [Acc24_m mean(Acc24)];
nmi24_m = [nmi24_m mean(nmi24)];
%SCML
%for N_iter = 1:5
     oldFolder = cd;
    cd ./SCML
    tx25 = tic;
    [label25] = sc_ml(M,m,lambda2);
    tEnd25 = toc(tx25);
    cd(oldFolder);
    Time25 = [Time25 tEnd25];
    if data=="aucs.mat" %if nargin==4
        label25(node_nolabels)=[]; %no consider nodes without label in calculation of accuracy
    end
    Acc25 = [Acc25 (n-wrong(c,label25'))/(n)];
    [T25] = confusion_matrix(c,label25');
    nmi25 = [nmi25 NMI(T25)];
%end
Time25_m = [Time25_m mean(Time25)];
Acc25_m = [Acc25_m mean(Acc25)];
nmi25_m = [nmi25_m mean(nmi25)];
%PMM
%for N_iter = 1:5
oldFolder = cd;
    cd ./PMM
    tx26 = tic;
    [label26] = PMM(A,  m, ell, maxKmeans);
    tEnd26 = toc(tx26);
    cd(oldFolder);
    Time26 = [Time26 tEnd26];
    if data=="aucs.mat" %if nargin==4
        label26(node_nolabels)=[]; %no consider nodes without label in calculation of accuracy
    end
    Acc26 = [Acc26 (n-wrong(c,label26'))/(n)];
    [T26] = confusion_matrix(c,label26');
    nmi26 = [nmi26 NMI(T26)];
%end
Time26_m = [Time26_m mean(Time26)];
Acc26_m = [Acc26_m mean(Acc26)];
nmi26_m = [nmi26_m mean(nmi26)];
end

GL=[mean(Time1_m);mean(Acc1_m);mean(nmi1_m)];
CVP1=[mean(Time2_m);mean(Acc2_m);mean(nmi2_m)];
CVP3=[mean(Time3_m);mean(Acc3_m);mean(nmi3_m)];
CVP5=[mean(Time4_m);mean(Acc4_m);mean(nmi4_m)];
CVP7=[mean(Time5_m);mean(Acc5_m);mean(nmi5_m)];
CVP9=[mean(Time6_m);mean(Acc6_m);mean(nmi6_m)];
MA_r2=[mean(Time7_m);mean(Acc7_m);mean(nmi7_m)];
MA_r3=[mean(Time8_m);mean(Acc8_m);mean(nmi8_m)];
MVP12=[mean(Time9_m);mean(Acc9_m);mean(nmi9_m)];
MVP13=[mean(Time10_m);mean(Acc10_m);mean(nmi10_m)];
MVP32=[mean(Time11_m);mean(Acc11_m);mean(nmi11_m)];
MVP33=[mean(Time12_m);mean(Acc12_m);mean(nmi12_m)];
MVP52=[mean(Time13_m);mean(Acc13_m);mean(nmi13_m)];
MVP53=[mean(Time14_m);mean(Acc14_m);mean(nmi14_m)];
MVP72=[mean(Time15_m);mean(Acc15_m);mean(nmi15_m)];
MVP73=[mean(Time16_m);mean(Acc16_m);mean(nmi16_m)];
MVP92=[mean(Time17_m);mean(Acc17_m);mean(nmi17_m)];
MVP93=[mean(Time18_m);mean(Acc18_m);mean(nmi18_m)];
CoRegm=[mean(Time19_m);mean(Acc19_m);mean(nmi19_m)];
AWPm=[mean(Time20_m);mean(Acc20_m);mean(nmi20_m)];
MCGCm=[mean(Time21_m);mean(Acc21_m);mean(nmi21_m)];
PM=[mean(Time22_m);mean(Acc22_m);mean(nmi22_m)];
MT=[mean(Time23_m);mean(Acc23_m);mean(nmi23_m)];
INFOMAP =[mean(Time24_m);mean(Acc24_m);mean(nmi24_m)];
SCML =[mean(Time25_m);mean(Acc25_m);mean(nmi25_m)];
PMMm =[mean(Time26_m);mean(Acc26_m);mean(nmi26_m)];

T=table(Methods,GL,CVP1,CVP3,CVP5,CVP7,CVP9,MA_r2,MA_r3,MVP12,MVP13,MVP32,MVP33,MVP52,MVP53,MVP72,MVP73,MVP92,MVP93,CoRegm,AWPm,MCGCm,PM,MT,INFOMAP,SCML,PMMm);
matrix_name = cellfun(@(s)s(1:end-4),cellstr(data),'uni',0);
file_name = ['Real_noise2_' matrix_name{1} '_' num2str(s) '.csv'];
writetable(T,file_name)
end