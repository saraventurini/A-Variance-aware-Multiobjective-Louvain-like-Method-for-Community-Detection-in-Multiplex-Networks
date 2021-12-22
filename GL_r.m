%Generalized Louvain 
%communities index randomly 
% Inputs : 
% M : 3 dim tensor, M(:,:,s)=adjacent matrix layer s  
% z : 1 = Recursive computation
%   : 0 = Just one level computation
%
% Output :
% COMTY, structure with the following information
% for each level i :
%   COMTY.COM{i} : vector of community partition
%   COMTY.SIZE{i} : vector of community sizes
%   COMTY.MOD{i} : vector of modularities of clustering on the layers 
%   COMTY.Average(i) : average of modularity on the layers 
%   COMTY.Niter(i) : Number of iteration before convergence
% 
function [COMTY, ending] = GL_r(M,z)
if nargin < 1
  error('not enough argument');
end

if nargin < 2
  z = 1;
end

S = size(M); 
N = S(1); %number of nodes
if length(S)==3
    k = S(3); %number of layers
else
    k = 1;
end

ending = 0;

%matrix without diagonal values
M2=zeros(N,N,k);
for s=1:k
    M2(:,:,s)=M(:,:,s)-diag(diag(M(:,:,s)));
end

%total weight of the graph 
m = zeros(k,1); 
for s=1:k
    if z==1
        M(:,:,s) = M(:,:,s) + diag(diag(M(:,:,s))); 
    end
    m(s) = sum(sum(M(:,:,s)))/2;
end
m2 = 2*m.^2;

Niter = 0; %number of iterations

if sum(m>0)==0 || N==1
  ending = 1;
  COMTY = 0;
  return;
end 
M(:,:,m==0) = [];
k = k - sum(m==0); 

COM = 1:S(1); % Community of node i

K(:,:) = sum(M(:,:,:),2); % K(i,s)= sum of wieght incident to node i in layer s 
SumIn = zeros(N,k);
Q = zeros(k,1);
for s=1:k
    SumIn(:,s) = diag(M(:,:,s)); %SumIn(i)= sum of weight inside community i (loops) 
                                  %At the beginning each node is a community
    Q(s) = compute_modularity(COM,M(:,:,s)); %modularity of layer s
end
SumTot = K; %SumIn(i)= sum of total weight of community i

Average = sum(Q)/k; %Average of modularity

NBc = cell(N,1);
%Neighbourhood 
for j=1:N
    for s=1:k
        NBc{j} = [NBc{j} find(M2(j,:,s))];
    end
    NBc{j} = unique(NBc{j});
end
%NBc{j}=neighbourhood of node j (sum of its neighbourhood on the layers)

gain = 1;
while (gain == 1)  %no increase of average possible moving one node 
  gain = 0;
  %Randomized labels (function reindex_com line 234)
  R = randperm(N);
  for g=1:N
    i = R(g);
    Ci = COM(i); %community of node i
    NB = NBc{i}; %neighbourhood of node i
    G = zeros(1,N); % Gain vector
    best_increase = -inf;
    Cnew = Ci;
    COM(i) = -1; %delate i from its neighbourhood
    %remove i from its community
    CNi = (COM==Ci); %list of nodes in community Ci, without i  
    tmp = sum(M(i,CNi,:),2); 
    GQ_i = (K(i,:).*SumTot(Ci,:))'./m2(:) - tmp(:)./m(:)- (K(i,:).^2)'./m2(:); %gain modularity due to remove i from its community
    SumTot(Ci,:) = SumTot(Ci,:) - K(i,:); %weights incident to Ci community 
    tmp2 = M(i,i,:);
    SumIn(Ci,:) = SumIn(Ci,:)' - 2*tmp(:) - tmp2(:); %sum of weight inside community Ci 
    
    for j=1:length(NB) %consider each neighbor j of node i
      Cj = COM(NB(j)); %community of of j
      if (G(Cj) == 0) %if we haven't already considereted the community of node j
        CNj = (COM==Cj); %nodes in community Cj, without j       
        tmp1 = sum(M(i,CNj,:),2);
        tmp2 = K(i,:);
        tmp3 = SumTot(Cj,:);
        DQ_vec = GQ_i(:) + tmp1(:)./m(:) - (tmp2(:).*tmp3(:))./m2(:); %modularity gain if node i is inserted into community Cj
        
        G(Cj) = sum(DQ_vec)/k; %Average gain 
        if G(Cj) > best_increase %if positive gain (average increases)
          best_increase = G(Cj); %gain average
          Q_t = DQ_vec; %gain of modularity;
          Cnew_t = Cj; %new community of node i
        end
      end
    end
    if best_increase > -10^(-15) %if best increase is positive
      Cnew = Cnew_t; %new community of node i 
      Q = Q + Q_t; %modularity
      Average = Average + best_increase; %average
    end
    %Recalculate
    Ck = (COM==Cnew);
    for s=1:k
        SumIn(Cnew,s) = SumIn(Cnew,s) + 2*sum(M(i,Ck,s)) + M(i,i,s);
        SumTot(Cnew,s) = SumTot(Cnew,s) + K(i,s);
    end
    COM(i) = Cnew; %insert node i in the new community 
    if (Cnew ~= Ci) %no increase of average possible moving one node 
      gain = 1;
    end
  end
  Niter = Niter + 1;
end
Niter = Niter - 1;
[COM] = reindex_com(COM); %reindex the communities
%output
COMTY.COM{1} = COM;
COMTY.MOD{1} = Q'; 
COMTY.Average(1) = Average;
COMTY.Niter(1) = Niter;


% Perform part 2
if (z == 1) 
  %matrix of the reduced network
  Mnew = M;  
  COMcur = COM; %current communities
  COMfull = COM; %communities in the original graph
  j = 2; %number of pass (step1+step2)
  while 1
    Mold = Mnew;
    S2 = size(Mold);
    Nnode = S2(1);
    
    COMu = unique(COMcur);
    Ncom = length(COMu);
    ind_com = sparse(Ncom,Nnode);
    ind_com_full = sparse(Ncom,N);
    for p=1:Ncom
      ind = find(COMcur==p);
      ind_com(p,1:length(ind)) = ind;
    end
    for p=1:Ncom
      ind = find(COMfull==p);
      ind_com_full(p,1:length(ind)) = ind;
    end
    Mnew = [];
    for s=1:k
        Mnew(:,:,s) = zeros(Ncom,Ncom); %new matrix (each node  is a community)
        for mm=1:Ncom
            for n=mm:Ncom
                ind1 = ind_com(mm,:);
                ind2 = ind_com(n,:);
                %weights of edges between communities
                 Mnew(mm,n,s) = sum(sum(Mold(ind1(ind1>0),ind2(ind2>0),s)));
                 Mnew(n,mm,s) = sum(sum(Mold(ind1(ind1>0),ind2(ind2>0),s)));
            end
        end
    end
    [COMt, e] = GL_r(Mnew,0); %apply the first phase on this reduced network 
    if (e ~= 1)
      COMfull = sparse(1,N);
      COMcur = COMt.COM{1};
      for p=1:Ncom
        ind1 = ind_com_full(p,:);
        COMfull(ind1(ind1>0)) = COMcur(p);
      end
      [COMfull2] = reindex_com(COMfull); %reindex the communitites
      %output
      COMTY.COM{j} = COMfull2;
      COMTY.MOD{j} = COMt.MOD{1};  
      COMTY.Average(j) = COMt.Average(1);
      COMTY.Niter(j) = COMt.Niter;
      Ind = (COMfull2 == COMTY.COM{j-1});
      if (sum(Ind) == length(Ind)) %no changes
        return;
      end
    else
      return;
    end
    j = j + 1;
  end
end
end

%Compute modulartiy
function MOD = compute_modularity(C,Mat)

m = sum(sum(Mat))/2; %total weight 
MOD = 0;
COMu = unique(C);

for j=1:length(COMu)
    Cj = (C==COMu(j)); %list of nodes in Cj
    Ec = sum(sum(Mat(Cj,Cj))); %sum of weights between nodes in Cj  
    Et = sum(sum(Mat(Cj,:))); %sum of weights of nodes incident in nodes of Cj 
    if Et>0
        MOD = MOD + Ec/(2*m)-(Et/(2*m))^2;
    end
end
end

%Re-index community partition randomly
function [C] = reindex_com(COMold)
[~,~,C]=unique(COMold);
end