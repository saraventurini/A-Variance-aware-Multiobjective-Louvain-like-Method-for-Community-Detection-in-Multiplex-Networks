% From adjecency matrix to file for INFOMAP method 
% First row: *Multilayer
% Each row: layer_id node_id layer_id node_id 
% INPUT: M adjacency matrix tensor. A(:,:,l) is adjacency matrix l-layer.
%        n number of nodes. 
%        k number of layers
function INFOMAP_adjacency_matrix(M,n,k)
oldFolder = cd;
fid = fopen('.\INFOMAP\adjacency_INFOMAP.net','w');
fprintf(fid,'*Multilayer \n');
%intraedges
for s=1:k %each layer
    for i=1:n %each node 
        for j=i+1:n
            if M(i,j)==1
                fprintf(fid,'%d %d %d %d \n',s,i,s,j);
                fprintf(fid,'%d %d %d %d \n',s,j,s,i); %symmetric matrix
            end
        end
    end
end

%interedges (for us: connecting same nodes on consecutive layers) (otherwise different communities on different layers)
for s=1:(k-1)
    for i=1:n
        fprintf(fid,'%d %d %d %d \n',s,i,s+1,i);
    end
end

fclose(fid);
cd(oldFolder)
end