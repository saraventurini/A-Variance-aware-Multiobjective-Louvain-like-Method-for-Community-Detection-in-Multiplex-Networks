%Tests Informative Case
%The function tests each methods on a SBM graph with 4 communities of size 125, p=0.1 probability arc intercommunity 
%and prints the average results on 10 runs in a csv file called "data.csv"
%INPUT: q probability arc intracommunity 
%       k number of layers

function run(q,k)

Methods={'Time';'Acc';'nmi'};

p=0.1; %probability arc intercommunity 

m=4; %number of community
L=125; %dimension of each community (every community same size)
n=m*L; %each layer size 
%Expected Communities
c=[];
for i=1:m
    for d=1:L
        c = [c i];
    end
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
addpath(genpath('ThePowermeanLaplacianForMultilayerGraphClustering'))
z = -10;

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

for N_iter2 = 1:10 %number of runs
%MAtrix 
M=adjacent_matrix_generator_multi(m,L,p,q,k);
for i=1:k
    A{i}=M(:,:,i);
end
%Run the methods on the MAtrix
    %Louvain Expansion Method Average
            tx1 = tic;
    [COMTY1, ending1] = EA_s(M); 
    tEnd1 = toc(tx1);
    if ending1==0
    Time1 = [Time1 tEnd1]; %time
    Acc1 = [Acc1 (n-wrong(c,COMTY1.COM{end}))/(n)]; %accuracy 
    [T1] = confusion_matrix(c,COMTY1.COM{end}'); %confusion MAtrix
    nmi1 = [nmi1 NMI(T1)]; %NMI
    end
    %Louvain Expansion Method Functin F-, lambda=0.1
            tx2 = tic;
    [COMTY2, ending2] = EVM(M,0.1);
    tEnd2 = toc(tx2);
    if ending2==0
    Time2 = [Time2 tEnd2];
    Acc2 = [Acc2 (n-wrong(c,COMTY2.COM{end}))/(n)];
    [T2] = confusion_matrix(c,COMTY2.COM{end}');
    nmi2 = [nmi2 NMI(T2)];
    end
    %Louvain Expansion Method Functin F-, lambda=0.3
        tx3 = tic;
    [COMTY3, ending3] = EVM(M,0.3); 
    tEnd3 = toc(tx3);
    if ending3==0
    Time3 = [Time3 tEnd3];
    Acc3 = [Acc3 (n-wrong(c,COMTY3.COM{end}))/(n)];
    [T3] = confusion_matrix(c,COMTY3.COM{end}');
    nmi3 = [nmi3 NMI(T3)];
    end
    %Louvain Expansion Method Functin F-, lambda=0.5
        tx4 = tic;
    [COMTY4, ending4] = EVM(M,0.5);  
    tEnd4 = toc(tx4);
    if ending4==0
    Time4 = [Time4 tEnd4];
    Acc4 = [Acc4 (n-wrong(c,COMTY4.COM{end}))/(n)];
    [T4] = confusion_matrix(c,COMTY4.COM{end}');
    nmi4 = [nmi4 NMI(T4)];
    end
    %Louvain Expansion Method Functin F-, lambda=0.7
        tx5 = tic;
    [COMTY5, ending5] = EVM(M,0.7);  
    tEnd5 = toc(tx5);
    if ending5==0
    Time5 = [Time5 tEnd5];
    Acc5 = [Acc5 (n-wrong(c,COMTY5.COM{end}))/(n)];
    [T5] = confusion_matrix(c,COMTY5.COM{end}');
    nmi5 = [nmi5 NMI(T5)];
    end
    %Louvain Expansion Method Functin F-, lambda=0.9
        tx6 = tic;
    [COMTY6, ending6] = EVM(M,0.9);  
    tEnd6 = toc(tx6);
    if ending6==0
    Time6 = [Time6 tEnd6];
    Acc6 = [Acc6 (n-wrong(c,COMTY6.COM{end}))/(n)];
    [T6] = confusion_matrix(c,COMTY6.COM{end}');
    nmi6 = [nmi6 NMI(T6)];
    end
    %Louvain Multiobjective Method Average, h=2
        tx7 = tic;
    [L7, ending7] = MA_s(M,2); 
    tEnd7 = toc(tx7);
    if ending7
    Time7 = [Time7 tEnd7];
    Acc7 = [Acc7 (n-wrong(c,L7{1}{2}{2}'))/(n)];
    [T7] = confusion_matrix(c,L7{1}{2}{2});
    nmi7 = [nmi7 NMI(T7)];
    end
    %Louvain Multiobjective Method Average, h=3
    tx8 = tic;
    [L8, ending8] = MA_s(M,3); 
    tEnd8 = toc(tx8);
    if ending8
    Time8 = [Time8 tEnd8];
    Acc8 = [Acc8 (n-wrong(c,L8{1}{2}{2}'))/(n)];
    [T8] = confusion_matrix(c,L8{1}{2}{2});
    nmi8 = [nmi8 NMI(T8)];
    end
    %Louvain Multiobjective Method Function F-, lambda=0.1, h=2
    tx9 = tic;
    [L9, ending9] = MVM(M,0.1,2); 
    tEnd9 = toc(tx9);
    if ending9
    Time9 = [Time9 tEnd9];
    Acc9 = [Acc9 (n-wrong(c,L9{1}{2}{2}'))/(n)];
    [T9] = confusion_matrix(c,L9{1}{2}{2});
    nmi9 = [nmi9 NMI(T9)];   
    end
    %Louvain Multiobjective Method Function F-, lambda=0.1, h=3
    tx10 = tic;
    [L10, ending10] = MVM(M,0.1,3); 
    tEnd10 = toc(tx10);
    if ending10
    Time10 = [Time10 tEnd10];
    Acc10 = [Acc10 (n-wrong(c,L10{1}{2}{2}'))/(n)];
    [T10] = confusion_matrix(c,L10{1}{2}{2});
    nmi10 = [nmi10 NMI(T10)];   
    end
    %Louvain Multiobjective Method Function F-, lambda=0.3, h=2
    tx11 = tic;
    [L11, ending11] = MVM(M,0.3,2); 
    tEnd11 = toc(tx11);
    if ending11
    Time11 = [Time11 tEnd11];
    Acc11 = [Acc11 (n-wrong(c,L11{1}{2}{2}'))/(n)];
    [T11] = confusion_matrix(c,L11{1}{2}{2});
    nmi11 = [nmi11 NMI(T11)];  
    end
    %Louvain Multiobjective Method Function F-, lambda=0.3, h=3
    tx12 = tic;
    [L12, ending12] = MVM(M,0.3,3); 
    tEnd12 = toc(tx12);
    if ending12
    Time12 = [Time12 tEnd12];
    Acc12 = [Acc12 (n-wrong(c,L12{1}{2}{2}'))/(n)];
    [T12] = confusion_matrix(c,L12{1}{2}{2});
    nmi12 = [nmi12 NMI(T12)]; 
    end
    %Louvain Multiobjective Method Function F-, lambda=0.5, h=2
    tx13 = tic;
    [L13, ending13] = MVM(M,0.5,2); 
    tEnd13 = toc(tx13);
    if ending13
    Time13 = [Time13 tEnd13];
    Acc13 = [Acc13 (n-wrong(c,L13{1}{2}{2}'))/(n)];
    [T13] = confusion_matrix(c,L13{1}{2}{2});
    nmi13 = [nmi13 NMI(T13)]; 
    end
    %Louvain Multiobjective Method Function F-, lambda=0.5, h=3
    tx14 = tic;
    [L14, ending14] = MVM(M,0.5,3); 
    tEnd14 = toc(tx14);
    if ending14
    Time14 = [Time14 tEnd14];
    Acc14 = [Acc14 (n-wrong(c,L14{1}{2}{2}'))/(n)];
    [T14] = confusion_matrix(c,L14{1}{2}{2});
    nmi14 = [nmi14 NMI(T14)];   
    end
    %Louvain Multiobjective Method Function F-, lambda=0.7, h=2
    tx15 = tic;
    [L15, ending15] = MVM(M,0.7,2); 
    tEnd15 = toc(tx15);
    if ending15
    Time15 = [Time15 tEnd15];
    Acc15 = [Acc15 (n-wrong(c,L15{1}{2}{2}'))/(n)];
    [T15] = confusion_matrix(c,L15{1}{2}{2});
    nmi15 = [nmi15 NMI(T15)];
    end
    %Louvain Multiobjective Method Function F-, lambda=0.7, h=3
    tx16 = tic;
    [L16, ending16] = MVM(M,0.7,3); 
    tEnd16 = toc(tx16);
    if ending16
    Time16 = [Time16 tEnd16];
    Acc16 = [Acc16 (n-wrong(c,L16{1}{2}{2}'))/(n)];
    [T16] = confusion_matrix(c,L16{1}{2}{2});
    nmi16 = [nmi16 NMI(T16)]; 
    end
    %Louvain Multiobjective Method Function F-, lambda=0.9, h=2
        tx17 = tic;
    [L17, ending17] = MVM(M,0.9,2); 
    tEnd17 = toc(tx17);
    if ending17
    Time17 = [Time17 tEnd17];
    Acc17 = [Acc17 (n-wrong(c,L17{1}{2}{2}'))/(n)];
    [T17] = confusion_matrix(c,L17{1}{2}{2});
    nmi17 = [nmi17 NMI(T17)];
    end
    %Louvain Multiobjective Method Function F-, lambda=0.9, h=3
    tx18 = tic;
    [L18, ending18] = MVM(M,0.9,3); 
    tEnd18 = toc(tx18);
    if ending18
    Time18 = [Time18 tEnd18];
    Acc18 = [Acc18 (n-wrong(c,L18{1}{2}{2}'))/(n)];
    [T18] = confusion_matrix(c,L18{1}{2}{2});
    nmi18 = [nmi18 NMI(T18)];
    end
    %CoReg
    tx19 = tic;
    [label19] = spectral_centroid_multiview_onkernel(A, k, m, lambda);
    tEnd19 = toc(tx19);
    Time19 = [Time19 tEnd19];
    Acc19 = [Acc19 (n-wrong(c,label19'))/(n)];
    [T19] = confusion_matrix(c,label19);
    nmi19 = [nmi19 NMI(T19)];
    %AWP
       tx20 = tic;
    for t=1:k
        embedding20{t} = spectral_embedding(A{t}, m);
    end
    label20 = AWP(embedding20, opts);
    tEnd20 = toc(tx20);
    Time20 = [Time20 tEnd20];
    Acc20 = [Acc20 (n-wrong(c,label20'))/(n)];
    [T20] = confusion_matrix(c,label20);
    nmi20 = [nmi20 NMI(T20)];
    %MCGC
        tx21 = tic;
    [label21] = obj_MVCC(M,k,m,beta);
    tEnd21 = toc(tx21);
    Time21 = [Time21 tEnd21];
    Acc21 = [Acc21 (n-wrong(c,label21'))/(n)];
    [T21] = confusion_matrix(c,label21);
    nmi21 = [nmi21 NMI(T21)];
    tx22 = tic;
    %PM
    [label22] = clustering_multilayer_graphs_with_power_mean_laplacian(A, z, m);
    tEnd22 = toc(tx22);
    Time22 = [Time22 tEnd22];
    Acc22 = [Acc22 (n-wrong(c,label22'))/(n)];
    [T22] = confusion_matrix(c,label22');
    nmi22 = [nmi22 NMI(T22)];
end 

%average values
CA=[mean(Time1);mean(Acc1);mean(nmi1)];
CVM1=[mean(Time2);mean(Acc2);mean(nmi2)];
CVM3=[mean(Time3);mean(Acc3);mean(nmi3)];
CVM5=[mean(Time4);mean(Acc4);mean(nmi4)];
CVM7=[mean(Time5);mean(Acc5);mean(nmi5)];
CVM9=[mean(Time6);mean(Acc6);mean(nmi6)];
MA2=[mean(Time7);mean(Acc7);mean(nmi7)];
MA3=[mean(Time8);mean(Acc8);mean(nmi8)];
MVM12=[mean(Time9);mean(Acc9);mean(nmi9)];
MVM13=[mean(Time10);mean(Acc10);mean(nmi10)];
MVM32=[mean(Time11);mean(Acc11);mean(nmi11)];
MVM33=[mean(Time12);mean(Acc12);mean(nmi12)];
MVM52=[mean(Time13);mean(Acc13);mean(nmi13)];
MVM53=[mean(Time14);mean(Acc14);mean(nmi14)];
MVM72=[mean(Time15);mean(Acc15);mean(nmi15)];
MVM73=[mean(Time16);mean(Acc16);mean(nmi16)];
MVM92=[mean(Time17);mean(Acc17);mean(nmi17)];
MVM93=[mean(Time18);mean(Acc18);mean(nmi18)];
CoRegm=[mean(Time19);mean(Acc19);mean(nmi19)];
AWPm=[mean(Time20);mean(Acc20);mean(nmi20)];
MCGCm=[mean(Time21);mean(Acc21);mean(nmi21)];
PM=[mean(Time22);mean(Acc22);mean(nmi22)];

T=table(Methods,CA,CVM1,CVM3,CVM5,CVM7,CVM9,MA2,MA3,MVM12,MVM13,MVM32,MVM33,MVM52,MVM53,MVM72,MVM73,MVM92,MVM93,CoRegm,AWPm,MCGCm,PM);
writetable(T,'data.csv') %csv file
end