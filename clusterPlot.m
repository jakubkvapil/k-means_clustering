function clusterPlot = clusterPlot(m,G,c)
Gsize= size(G);
k=Gsize(1,2);
Gmodif=zeros(length(G(:,1)),length(G(1,:)));
for i = 1: length(G(:,1))
    for j = 1: length(G(1,:))
        Gmodif(i,j)=G(i,j)*j;
    end
end
idx=sum(Gmodif,2);
figure;
for i = 1:k
    scatter(m(idx==i,1),m(idx==i,2),'filled');
    hold on
end
scatter(c(:,1),c(:,2),200,'x', 'k');
title 'Cluster Assignments and Centroids';
xlabel('x1');
ylabel('x2');