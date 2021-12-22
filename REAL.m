%Real Datasets-Informative case
%The function tests each methods and prints the results in a csv file called "data.csv"
%INPUT: W_cell cell containing adjacent matrix of each layer
%       labels vector right communities 
%       data name dataset
%       node_nolabels (just for aucs dataset) list of nodes with no true community assigned, then no consider for accuracy calculation 
%Choose between lines 23 and 24 depending on the dataset 

function REAL(W_cell,labels,data,node_nolabels)

Methods={'Time';'Acc';'NMI'};

%Matrix
N = size(W_cell{1},1); %nodes
if data=="aucs.mat" %if nargin==4
    n = N-length(node_nolabels); %no consider nodes without label in calculation of accuracy
else
    n = N;
end
k = size(W_cell,1); %layers
for i=1:k
    M(:,:,i)=full(W_cell{i}); %Matrix
end
for i=1:k
    A{i}=M(:,:,i);
end
%file adjacency matrix for MultiTensor method
MULTITENSOR_adjacency_matrix(M,n);
%file adjacency matrix for INFOMAP method "adjacency_INFOMAP.net"
INFOMAP_adjacency_matrix(M,N,k); 

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

    %Generalized Louvain
     tx1 = tic;
    [COMTY1, ending1] = GL_s(M); 
    tEnd1 = toc(tx1);
    if ending1==0
    Time1 = tEnd1; %time
    if data=="aucs.mat" %if nargin==4
        COMTY1.COM{end}(node_nolabels)=[]; %no consider nodes without label in calculation of accuracy
    end 
    Acc1 = (n-wrong(c,COMTY1.COM{end}))/(n); %accuracy 
    [T1] = confusion_matrix(c,COMTY1.COM{end}'); %confusion Matrix
    nmi1 = NMI(T1); %NMI
    end
GL=[Time1;Acc1;nmi1];

    %Louvain Expansion Method Functin F-, lambda=0.1
            tx2 = tic;
    [COMTY2, ending2] = EVM(M,0.1);
    tEnd2 = toc(tx2);
    if ending2==0
    Time2 = tEnd2;
    if data=="aucs.mat" %if nargin==4
        COMTY2.COM{end}(node_nolabels)=[]; %no consider nodes without label in calculation of accuracy
    end 
    Acc2 = (n-wrong(c,COMTY2.COM{end}))/(n);
    [T2] = confusion_matrix(c,COMTY2.COM{end}');
    nmi2 = NMI(T2);
    end
CVM1=[Time2;Acc2;nmi2];

    %Louvain Expansion Method Functin F-, lambda=0.3
        tx3 = tic;
    [COMTY3, ending3] = EVM(M,0.3); 
    tEnd3 = toc(tx3);
    if ending3==0
    Time3 = tEnd3;
    if data=="aucs.mat" %if nargin==4
        COMTY3.COM{end}(node_nolabels)=[]; %no consider nodes without label in calculation of accuracy
    end 
    Acc3 = (n-wrong(c,COMTY3.COM{end}))/(n);
    [T3] = confusion_matrix(c,COMTY3.COM{end}');
    nmi3 = NMI(T3);
    end
CVM3=[Time3;Acc3;nmi3];

    %Louvain Expansion Method Functin F-, lambda=0.5
        tx4 = tic;
    [COMTY4, ending4] = EVM(M,0.5);  
    tEnd4 = toc(tx4);
    if ending4==0
    Time4 = tEnd4;
    if data=="aucs.mat" %if nargin==4
        COMTY4.COM{end}(node_nolabels)=[]; %no consider nodes without label in calculation of accuracy
    end 
    Acc4 = (n-wrong(c,COMTY4.COM{end}))/(n);
    [T4] = confusion_matrix(c,COMTY4.COM{end}');
    nmi4 = NMI(T4);
    end
CVM5=[Time4;Acc4;nmi4];

    %Louvain Expansion Method Functin F-, lambda=0.7
        tx5 = tic;
    [COMTY5, ending5] = EVM(M,0.7);  
    tEnd5 = toc(tx5);
    if ending5==0
    Time5 = tEnd5;
    if data=="aucs.mat" %if nargin==4
        COMTY5.COM{end}(node_nolabels)=[]; %no consider nodes without label in calculation of accuracy
    end 
    Acc5 = (n-wrong(c,COMTY5.COM{end}))/(n);
    [T5] = confusion_matrix(c,COMTY5.COM{end}');
    nmi5 = NMI(T5);
    end
CVM7=[Time5;Acc5;nmi5];

    %Louvain Expansion Method Functin F-, lambda=0.9
        tx6 = tic;
    [COMTY6, ending6] = EVM(M,0.9);  
    tEnd6 = toc(tx6);
    if ending6==0
    Time6 = tEnd6;
    if data=="aucs.mat" %if nargin==4
        COMTY6.COM{end}(node_nolabels)=[]; %no consider nodes without label in calculation of accuracy
    end 
    Acc6 = (n-wrong(c,COMTY6.COM{end}))/(n);
    [T6] = confusion_matrix(c,COMTY6.COM{end}');
    nmi6 = NMI(T6);
    end
CVM9=[Time6;Acc6;nmi6];

   %Louvain Multiobjective Method Average, h=2
        tx7 = tic;
    [L7, ending7] = MA_s(M,2); 
    tEnd7 = toc(tx7);
    if ending7
    Time7 = tEnd7;
    if data=="aucs.mat" %if nargin==4
        L7{1}{2}{2}(node_nolabels)=[]; %no consider nodes without label in calculation of accuracy
    end 
    Acc7 = (n-wrong(c,L7{1}{2}{2}'))/(n);
    [T7] = confusion_matrix(c,L7{1}{2}{2});
    nmi7 = NMI(T7);
    end
MA_s2=[Time7;Acc7;nmi7];
    %Louvain Multiobjective Method Average, h=3
    tx8 = tic;
    [L8, ending8] = MA_s(M,3); 
    tEnd8 = toc(tx8);
    if ending8
    Time8 = tEnd8;
    if data=="aucs.mat" %if nargin==4
        L8{1}{2}{2}(node_nolabels)=[]; %no consider nodes without label in calculation of accuracy
    end 
    Acc8 = (n-wrong(c,L8{1}{2}{2}'))/(n);
    [T8] = confusion_matrix(c,L8{1}{2}{2});
    nmi8 = NMI(T8);
    end
MA_s3=[Time8;Acc8;nmi8];

    %Louvain Multiobjective Method Function F-, lambda=0.1, h=2
    tx9 = tic;
    [L9, ending9] = MVM(M,0.1,2); 
    tEnd9 = toc(tx9);
    if ending9
    Time9 = tEnd9;
    if data=="aucs.mat" %if nargin==4
        L9{1}{2}{2}(node_nolabels)=[]; %no consider nodes without label in calculation of accuracy
    end
    Acc9 = (n-wrong(c,L9{1}{2}{2}'))/(n);
    [T9] = confusion_matrix(c,L9{1}{2}{2});
    nmi9 = NMI(T9);   
    end
MVM12=[Time9;Acc9;nmi9];
    keyboard
    %Louvain Multiobjective Method Function F-, lambda=0.1, h=3
    tx10 = tic;
    [L10, ending10] = MVM(M,0.1,3); 
    tEnd10 = toc(tx10);
    if ending10
    Time10 = tEnd10;
    if data=="aucs.mat" %if nargin==4
        L10{1}{2}{2}(node_nolabels)=[]; %no consider nodes without label in calculation of accuracy
    end
    Acc10 = (n-wrong(c,L10{1}{2}{2}'))/(n);
    [T10] = confusion_matrix(c,L10{1}{2}{2});
    nmi10 = NMI(T10);   
    end
MVM13=[Time10;Acc10;nmi10];

    %Louvain Multiobjective Method Function F-, lambda=0.3, h=2
    tx11 = tic;
    [L11, ending11] = MVM(M,0.3,2); 
    tEnd11 = toc(tx11);
    if ending11
    Time11 = tEnd11;
    if data=="aucs.mat" %if nargin==4
        L11{1}{2}{2}(node_nolabels)=[]; %no consider nodes without label in calculation of accuracy
    end
    Acc11 = (n-wrong(c,L11{1}{2}{2}'))/(n);
    [T11] = confusion_matrix(c,L11{1}{2}{2});
    nmi11 = NMI(T11);  
    end
MVM32=[Time11;Acc11;nmi11];

    %Louvain Multiobjective Method Function F-, lambda=0.3, h=3
    tx12 = tic;
    [L12, ending12] = MVM(M,0.3,3); 
    tEnd12 = toc(tx12);
    if ending12
    Time12 = tEnd12;
    if data=="aucs.mat" %if nargin==4
        L12{1}{2}{2}(node_nolabels)=[]; %no consider nodes without label in calculation of accuracy
    end
    Acc12 = (n-wrong(c,L12{1}{2}{2}'))/(n);
    [T12] = confusion_matrix(c,L12{1}{2}{2});
    nmi12 = NMI(T12); 
    end
MVM33=[Time12;Acc12;nmi12];

    %Louvain Multiobjective Method Function F-, lambda=0.5, h=2
    tx13 = tic;
    [L13, ending13] = MVM(M,0.5,2); 
    tEnd13 = toc(tx13);
    if ending13
    Time13 = tEnd13;
    if data=="aucs.mat" %if nargin==4
        L13{1}{2}{2}(node_nolabels)=[]; %no consider nodes without label in calculation of accuracy
    end
    Acc13 = (n-wrong(c,L13{1}{2}{2}'))/(n);
    [T13] = confusion_matrix(c,L13{1}{2}{2});
    nmi13 = NMI(T13); 
    end
MVM52=[Time13;Acc13;nmi13];

    %Louvain Multiobjective Method Function F-, lambda=0.5, h=3
    tx14 = tic;
    [L14, ending14] = MVM(M,0.5,3); 
    tEnd14 = toc(tx14);
    if ending14
    Time14 = tEnd14;
    if data=="aucs.mat" %if nargin==4
        L14{1}{2}{2}(node_nolabels)=[]; %no consider nodes without label in calculation of accuracy
    end
    Acc14 = (n-wrong(c,L14{1}{2}{2}'))/(n);
    [T14] = confusion_matrix(c,L14{1}{2}{2});
    nmi14 = NMI(T14);   
    end
MVM53=[Time14;Acc14;nmi14];

    %Louvain Multiobjective Method Function F-, lambda=0.7, h=2
    tx15 = tic;
    [L15, ending15] = MVM(M,0.7,2); 
    tEnd15 = toc(tx15);
    if ending15
    Time15 = tEnd15;
    if data=="aucs.mat" %if nargin==4
        L15{1}{2}{2}(node_nolabels)=[]; %no consider nodes without label in calculation of accuracy
    end
    Acc15 = (n-wrong(c,L15{1}{2}{2}'))/(n);
    [T15] = confusion_matrix(c,L15{1}{2}{2});
    nmi15 = NMI(T15);
    end
MVM72=[Time15;Acc15;nmi15];

    %Louvain Multiobjective Method Function F-, lambda=0.7, h=3
    tx16 = tic;
    [L16, ending16] = MVM(M,0.7,3); 
    tEnd16 = toc(tx16);
    if ending16
    Time16 = tEnd16;
    if data=="aucs.mat" %if nargin==4
        L16{1}{2}{2}(node_nolabels)=[]; %no consider nodes without label in calculation of accuracy
    end
    Acc16 = (n-wrong(c,L16{1}{2}{2}'))/(n);
    [T16] = confusion_matrix(c,L16{1}{2}{2});
    nmi16 = NMI(T16); 
    end
MVM73=[Time16;Acc16;nmi16];

    %Louvain Multiobjective Method Function F-, lambda=0.9, h=2
        tx17 = tic;
    [L17, ending17] = MVM(M,0.9,2); 
    tEnd17 = toc(tx17);
    if ending17
    Time17 = tEnd17;
    if data=="aucs.mat" %if nargin==4
        L17{1}{2}{2}(node_nolabels)=[]; %no consider nodes without label in calculation of accuracy
    end
    Acc17 = (n-wrong(c,L17{1}{2}{2}'))/(n);
    [T17] = confusion_matrix(c,L17{1}{2}{2});
    nmi17 = NMI(T17);
    end
MVM92=[Time17;Acc17;nmi17];

    %Louvain Multiobjective Method Function F-, lambda=0.9, h=3
    tx18 = tic;
    [L18, ending18] = MVM(M,0.9,3); 
    tEnd18 = toc(tx18);
    if ending18
    Time18 = tEnd18;
    if data=="aucs.mat" %if nargin==4
        L18{1}{2}{2}(node_nolabels)=[]; %no consider nodes without label in calculation of accuracy
    end
    Acc18 = (n-wrong(c,L18{1}{2}{2}'))/(n);
    [T18] = confusion_matrix(c,L18{1}{2}{2});
    nmi18 = NMI(T18);
    end
MVM93=[Time18;Acc18;nmi18];

    %CoReg
    tx19 = tic;
    [label19] = spectral_centroid_multiview_onkernel(A, k, m, lambda);
    tEnd19 = toc(tx19);
    Time19 = tEnd19;
    if data=="aucs.mat" %if nargin==4
        label19(node_nolabels)=[]; %no consider nodes without label in calculation of accuracy
    end
    Acc19 = (n-wrong(c,label19'))/(n);
    [T19] = confusion_matrix(c,label19);
    nmi19 = NMI(T19);
CoRegm=[Time19;Acc19;nmi19];

    %AWP
       tx20 = tic;
    for t=1:k
        embedding20{t} = spectral_embedding(A{t}, m);
    end
    label20 = AWP(embedding20, opts);
    tEnd20 = toc(tx20);
    Time20 = tEnd20;
    if data=="aucs.mat" %if nargin==4
        label20(node_nolabels)=[]; %no consider nodes without label in calculation of accuracy
    end
    Acc20 = (n-wrong(c,label20'))/(n);
    [T20] = confusion_matrix(c,label20);
    nmi20 = NMI(T20);
AWPm=[Time20;Acc20;nmi20];

%MCGC
        tx21 = tic;
    [label21] = obj_MVCC(M,k,m,beta);
    tEnd21 = toc(tx21);
    Time21 = tEnd21;
    if data=="aucs.mat" %if nargin==4
        label21(node_nolabels)=[]; %no consider nodes without label in calculation of accuracy
    end
    Acc21 = (n-wrong(c,label21'))/(n);
    [T21] = confusion_matrix(c,label21);
    nmi21 = NMI(T21);
    tx22 = tic;
MCGCm=[Time21;Acc21;nmi21];

    %PM
    [label22] = clustering_multilayer_graphs_with_power_mean_laplacian(A, z, m);
    tEnd22 = toc(tx22);
    Time22 = tEnd22;
    if data=="aucs.mat" %if nargin==4
        label22(node_nolabels)=[]; %no consider nodes without label in calculation of accuracy
    end
    Acc22 = (n-wrong(c,label22'))/(n);
    [T22] = confusion_matrix(c,label22');
    nmi22 = NMI(T22);
PM=[Time22;Acc22;nmi22];

    %MT
    %create file adjacency matrix
    MULTITENSOR_adjacency_matrix(M,N)
    oldFolder = cd;
    cd ./MultiTensor-master/python 
    %run the method   
    python_file_name = ['python main.py -k=' num2str(m) ' -l=' num2str(k) ' -a="adjacency.dat" -E="_endfile.dat"'];    
    tx23 = tic;
    system(python_file_name);
    Time23 = toc(tx23);
    %read the result
    cd(oldFolder)
    file_results = ['./MultiTensor-master/data/u_K' num2str(m) '_endfile.dat'];
    B = readtable(file_results);
    B = B{:,:};
    B(:,1)=[];
    [~, label23] = max(B,[],2);
    if data=="aucs.mat" %if nargin==4
        label23(node_nolabels)=[]; %no consider nodes without label in calculation of accuracy
    end
    Acc23 = (n-wrong(c,label23'))/(n);
    [T23] = confusion_matrix(c,label23');
    nmi23 = NMI(T23);
MT=[Time23;Acc23;nmi23];
    
    %INFOMAP
    %run the method
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
INFOMAP=[Time24;Acc24;nmi24];

    %SCML
    %run the method
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
SCML =[Time25;Acc25;nmi25];

    %PMM
    %run the method
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
PMMm =[Time26;Acc26;nmi26];

T=table(Methods,GL,CVM1,CVM3,CVM5,CVM7,CVM9,MA_s2,MA_s3,MVM12,MVM13,MVM32,MVM33,MVM52,MVM53,MVM72,MVM73,MVM92,MVM93,CoRegm,AWPm,MCGCm,PM,MT,INFOMAP,SCML,PMMm);
matrix_name = cellfun(@(s)s(1:end-4),cellstr(data),'uni',0);
file_name = ['Real_info_' matrix_name{1} '.csv'];
writetable(T,file_name)

end