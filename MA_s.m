%Louvain Multiobjective Method Average 
%communities index by size
% Inputs :  
% M : 3 dim tensor, M(:,:,s)=adjacent matrix layer s  
% h : length of the filter 
% z : 1 = Recursive computation
%   : 0 = Just one level computation
%
% Output :
% L cell with the following information
%   L{1}=Q cell with modularity of each layer  
%   L{2}={COMcur COMfull COMindex} communities in the current graph, in  the
%   original graph, communities in the current graph reindexed 
%   L{3}=SumTot (sum of weights of the links incident to nodes in a community)
%   L{4}=SumIn  (sum of weights of the links inside a community)
%   L{5}=Average  (average)
%

function [L, ending] = MA_s(M,h,z)

if nargin < 1
  error('not enough argument');
end
if nargin < 2
  error('not h defined');
end
if nargin < 3
  z = 1;
end

S = size(M); 
N = S(1); %number of nodes
if length(S)==3
    k = S(3); %number of layers
else k = 1;
end 

ending = 1;

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

Niter = 0; %number of iterations

% Calculation of the beginning values
COM{1} = 1:N; %current community %At the beginning each node is a community 
COM{2} = 1:N; %original graph commuity 
COM{3} = 1:N; %community reindexed
%COM(i)= Community of node i

K(:,:) = sum(M(:,:,:),2); % K(i,s)= sum of wieght incident to node i in layer s  
SumIn = zeros(N,k);
Q = zeros(k,1);
for s=1:k   
    SumIn(:,s) = diag(M(:,:,s)); %SumIn(i)= sum of weight inside community i (loops) 
                                  %At the beginning each node is a community
    Q(s) = compute_modularity(COM{1},M(:,:,s)); %modularity of layer s
end
SumTot = K; %SumIn(i)= sum of total weight of community i

%average
Average = sum(Q)/k;

%If no edges in any layer or just one node
if sum(m>0)==0 || N==1
  ending = 0;
  L={[Q,{COM},SumTot,SumIn,Average]};
else
%Delate layers with m=0 
M(:,:,m==0) = [];
k = k - sum(m==0);
Q(m==0)=[];
m(m==0)=[];

%Initialization filter
L={[Q,{COM},SumTot,SumIn,Average]};

NBc = cell(N,1);
%Neighbourhood 
for j=1:N
    for s=1:k
        NBc{j} = [NBc{j} find(M2(j,:,s))];
    end
    NBc{j} = unique(NBc{j});
end
%NBc{j}=neighbourhood of node j (sum of its neighbourhood on the layers)


GAIN = 1;
while (GAIN == 1)  %filter has changed  
  L_old=L; %memorize initial filter 
  L_new=L; %new filter that are going to change during the phase
  %Index communitites by size (function reindex_com line 387)
  for i=1:N %for each node, in this order
      L_o=L; %memorize the initial filter 
      for l=1:length(L_o) %for each partition in the filter
          %delate the point from the filter
          L_new_o=L_new;
          index = cellfun(@(x) isequal(x,L_o{l}), L_new, 'UniformOutput', 1);
          L_new(index)=[];
          %apply step1 to this partition
          [L_new,U] = step_1 (L_o{l}{1},L_o{l}(2),L_o{l}{3},L_o{l}{4},L_o{l}{5},k,NBc,K,M,N,m,L_new,h,i);
          %if no new point, insert again the initial point in the filter
          if U==0
            L_new=L_new_o;
          end
      end
      L_new = cut_filter(L_new,h); %cut the filter at length h 
      L=L_new; %new filter
  end
  %Check if nothing happened
  if length(L_old)~=length(L_new)
      GAIN=1;
  else
      g = length(L_old);
      C_old=zeros(g,k);
      C_new=zeros(g,k);
      for l=1:length(L_old)
          C_old(l,:) = L_old{l}{1}';
          C_new(l,:) = L_new{l}{1}';
      end
      C = ismembertol(C_old,C_new,10^(-3));
      if size(C,1)*size(C,2) == sum(sum(C))
          GAIN=0; %no filter changes  
      end 
  end
  Niter = Niter + 1;
end
L = cut_filter(L,1); %one final element
L{1}{2}{3}=reindex_com(L{1}{2}{1})'; %reindex communities
end

% Perform part 2
if (z==1)
      %matrix of the reduced network
      Mnew=M; 
      COMcur = L{1}{2}{3}; %communities in the graph of this iteration
      COMfull = L{1}{2}{3}; %communities in the original graph
      LL=L; %memorize initial filter
      j=2;%number of pass (step1+step2) 
      S=1;
      while S %element of the filter has changed
        L_old=LL; %memorize filter
        Mold = Mnew;
        S2 = size(Mold); %number of nodes 
        Nnode = S2(1);
        
        COMu = unique(COMcur); %list of communities without repetitions and in sorted order
        Ncom = length(COMu); %number of communities
        ind_com = sparse(Ncom,Nnode); %zero matrix with a row for each community and a column for each node in this configuration
        ind_com_full =sparse(Ncom,N); %zero matrix with a row for each community and a column for each node of the original graph
    
        for p=1:Ncom %for each community
            ind = find(COMcur==p);
            ind_com(p,1:length(ind)) = ind; %in row p put the indeces of the nodes in community p in this configuration
        end
        for p=1:Ncom
            ind = find(COMfull==p);
            ind_com_full(p,1:length(ind)) = ind; %in row p put the indeces of the nodes in community p in the original graph 
        end
        Mnew=[];
        for s=1:k
            Mnew(:,:,s) = sparse(Ncom,Ncom); %new matrix (each node  is a community)
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
        %apply first step to this new matrix 
        [LL e] = MA_s(Mnew,h,0);
            COMfull = sparse(1,N); %communities of the original graph
            COMcur = LL{1}{2}{3}; %communities
            for p=1:Ncom 
                ind1 = ind_com_full(p,:); %nodes in community p in the original graph
                COMfull(ind1(ind1>0)) = COMcur(p); %community now of node i 
            end
            [COMfull] = reindex_com(COMfull); %reindex communitites of original graph
            LL{1}{2}{2}=COMfull;
            %communities do not change
            if isequal(L_old{1}{2}{2},LL{1}{2}{2}) %filter no changes
                L=LL;
                S=0;
            end             
        j = j + 1; %start another pass
      end
end
end

%Compute step_1 
function [L,U] = step_1(Q,COM,SumTot,SumIn,Average,k,NBc,K,M,N,m,L,h,i)
U=0; %if U=1 add at least an element to the list, otherwise insert again the initial point in the filter
COMcur=COM{1}{1}; %current community 
Ci = COMcur(i); %community of node i 

m2 = 2*m.^2;
NB = NBc{i}; %neighbourhood of node i

%reallocate all the values becouse they change in this loop but they are
%used also for the other partition of the list 
SumTot2=SumTot;
SumIn2=SumIn;

%remove i from its community 
COMcur2 = COMcur;
COMcur2(i) = -1; %delate i from its neighbourhood
CNi = (COMcur2==Ci); %list of nodes in Ci community, without i  
tmp = sum(M(i,CNi,:),2); 
GQ_i = (K(i,:).*SumTot(Ci,:))'./m2(:) - tmp(:)./m(:)- (K(i,:).^2)'./m2(:); %gain modularity due to remove i from its community
SumTot2(Ci,:) = SumTot2(Ci,:) - K(i,:); %weights incident to Ci community 
tmp2 = M(i,i,:);
SumIn2(Ci,:) = SumIn2(Ci,:)' - 2*tmp(:) - tmp2(:); %sum of weight inside community Ci
    

G = sparse(1,N);
for j=1:length(NB) %consider each neighbor j of node i
    SumTot3=SumTot2; %reallocate values
    SumIn3=SumIn2;
    Cj = COMcur(NB(j)); %community of node j
    if (G(Cj) == 0) %If not tried with another node of community Cj yet
        G(Cj)=1;
        COMcur3 = COMcur;
        COMcur3(i) = -1;
        CNj = (COMcur3==Cj);
        %Insert i in the community of j for each layer
        tmp1 = sum(M(i,CNj,:),2);
        tmp2 = K(i,:);
        tmp3 = SumTot3(Cj,:);
        DQ_vec = GQ_i(:) + tmp1(:)./m(:) - (tmp2(:).*tmp3(:))./m2(:); %modularity gain if node i is inserted into community Cj
    %average
    GA_j = sum(DQ_vec)/k; %gain average
    if GA_j > -10^(-14) %if positive gain (average increases)
        %Recalculate
        A_j = Average + GA_j; %average
        SumTot3(Cj,:) = SumTot3(Cj,:) + K(i,:); %SumTot
        tmp1 = sum(M(i,CNj,:),2);
        tmp2 = M(i,i,:);
        SumIn3(Cj,:) = SumIn3(Cj,:) + 2.*tmp1(:)' + tmp2(:)'; %SumIn
        COMcur4=COMcur;
        COMcur4(i)=Cj; %insert node i in the new community 
        COM{1}{1}=COMcur4;
        Q_j = Q + DQ_vec; %modularity 
        %add the partition to the filter
        [L] = add_to_filter(Q_j,COM,SumTot3,SumIn3,A_j,L,k);
        U=1;      
    end
    end
end
end

%add a point to the filter 
function [L] = add_to_filter(QQ,COM,SumTot,SumIn,Average,L,k)

if  ~isempty(L)

%1. Check if the new point is dominated by a point in the filter
DD=0;
l=1;
while (DD==0) && l<=length(L) 
        D=0;
        s=1;
        while (D==0) && s<=k
                if QQ(s)>L{l}{1}(s)
                    D=1;
                end
                s=s+1;
        end
        if ~D %the new point is dominated by L{l} 
            DD=1;
        end
        l=l+1;
end
if  ~DD %the new point is not dominated  
  DDD=0;
  l=1;
  while (DDD==0) && l<=length(L)
      D=1;
      s=1;
      while D && s<=k
          if abs(L{l}{1}(s) - QQ(s)) > 10^(-14) 
              D=0;
          end
          s = s + 1;
      end
      if D
          DDD=1; %same point
      end
      l = l + 1;
  end
  
  if DDD==0 
      
  %2.Check if the other points in the filter are dominated by the new point 

  for l=1:length(L)
        D=1;
        s=1;
        while D && s<=k
                if L{l}{1}(s)>QQ(s)+10^(-14)
                    D=0;
                end
                s=s+1;
        end
        if D %L{l} is dominated by the new point so I remove it from the filter 
            L{l}=[];
        end
  end
 
  empties =(cellfun(@isempty,L));
  L(empties) = [];
  L{end+1}=[QQ,COM,SumTot,SumIn,Average];  %add the new point to the filter
   
end
end
else %empty filter
    L={[QQ,COM,SumTot,SumIn,Average]};  %add the new point to the filter
end
end


%cut the filter
function [L] = cut_filter(L,h)

%average of every element of the filter 
if length(L)>h   
   Average=cell(1,length(L));
   for l=1:length(L)
        Average{l} = L{l}{5};
   end  
   
%remove the situations in the filter with the lower values
for t=1:(length(L)-h)
  [mn]=min([Average{:}]);
  v=find([Average{:}]==mn); %if more element with same average values, delate the last one 
                            %(therefore with h=1 I have same result as Louvain Expansion Method Average)
  idx_t=v(end);
  Average{idx_t}=Inf;
  L{idx_t}=[];
end

empties = cellfun(@isempty,L);
L(empties) = [];
end   
end
 
%Compute modulartiy
function MOD = compute_modularity(C,Mat)

S = size(Mat);
N = S(1);

m = sum(sum(Mat))/2; %total weight 

MOD = 0;
COMu = unique(C); %list of communities without repetiotions and in sorted order 
%for each community calculate modularity and then sum all together 
for j=1:length(COMu)
    Cj = (C==COMu(j)); %list of nodes in Cj
    Ec = sum(sum(Mat(Cj,Cj))); %sum of weights between nodes in Cj 
    Et = sum(sum(Mat(Cj,:))); %sum of weights of nodes incident in nodes of Cj 
    if Et>0
        MOD = MOD + Ec/(2*m)-(Et/(2*m))^2;
    end
end

end

% Re-index community partition by size 
function [C] = reindex_com(COMold)
C = sparse(1,length(COMold));
COMolds=sort(COMold);
COMu=COMolds([true;diff(COMolds(:))>0]);
S = sparse(1,length(COMu));
for l=1:length(COMu)
    S(l) = length(COMold(COMold==COMu(l)));
end
[Ss INDs] = sort(S,'descend');
for l=1:length(COMu)
    C(COMold==COMu(INDs(l))) = l;
end
end
