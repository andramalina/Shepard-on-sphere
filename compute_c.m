function cc=compute_c(nodes)
N=length(nodes);
dist=[];
for i=1:N
dist_i = real(acos(nodes(i,:)*nodes'));
  [sorted_dist,I] = sort(dist_i);
dist=[dist, sorted_dist(2)];
end
d=sum(dist)/N;
cc=1/(0.815*d);
end