%====================================================
load('5000HaltonNodesb2b3.mat');
nodes=[nodes(:,1), nodes(:,2), nodes(:,3)];
% load('1000spiralpoints.mat');
% nodes=[points(:,1), points(:,2), points(:,3)];
%====================================================
load('600spiralpoints.mat');
x = [points(:,1), points(:,2), points(:,3)];
% load('1000points-random.mat')
% random_points=points;
% TH_points = 2*pi*random_points;
% PH_points = asin(-1+2*random_points);
% [xP,yP,zP] = sph2cart(TH_points,PH_points,1);
% x = [xP' yP' zP'];
%====================================================
tic
nr_nodes=length(nodes);
nr_points = length(points);
Nw=8;   
Nz=17;
p=2;
c=compute_c(nodes)
%====================================================
result_for_rmse=zeros(1,nr_points);
result_for_mae=zeros(1,nr_points);
true_vals=zeros(1,nr_points);
%====================================================
for j=1:nr_points
    a=0; b=0;
    a=shep_modif_sphere(x(j,:), nodes, nr_nodes, Nw, Nz,p,c);
    b=func(x(j,:));
    result_for_mae(j) = abs(a-b);  
    result_for_rmse(j)=a;
    true_vals(j)=b;
end
%====================================================
mae=max(result_for_mae)
rmse=(sum((result_for_rmse-true_vals).^2)/nr_points)^1/2
figure(1)
scatter3(x(:,1),x(:,2),x(:,3),20,true_vals', 'filled')
colorbar
figure(2)
scatter3(x(:,1),x(:,2),x(:,3),20,result_for_rmse', 'filled')
colorbar
toc
%====================================================
function s_sum = shep_modif_sphere(x, nodes, nr_nodes, Nw, Nz, p, c)
f = zeros(nr_nodes, 1);
for j=1:nr_nodes
    f(j) = compute_Zj(j,nodes, Nz, x,c);
end
s = zeros(1, nr_nodes);
s_sum=Shepard_modif_triv(x,nodes,f,Nw,p); 
end
%==================================================================
function z_j =compute_Zj(j, nodes, Nz, x, c)
node_j = nodes(j,:);
if x==node_j
    z_j = func(node_j);
    return
end
neighb_j = find_neighbours(node_j, nodes, Nz);
a_j = find_coeff(neighb_j, Nz,c );
z_j =0;
for i=1:Nz
    t = real(geodesic_dist(x, neighb_j(i,:)));
    z_j = z_j + a_j(i) * phi(t,c);
end
z_j=z_j + a_j(Nz+1).*x(1) + a_j(Nz+2).*x(2) + a_j(Nz+3).*x(3);
end
%==================================================================
function a = find_coeff(nodes, Nz, c)
a=zeros(Nz+3,1);
M=zeros(Nz+3); 
f=zeros(Nz+3,1);
for i=1:Nz
    for j=1:Nz
            t = (geodesic_dist(nodes(i,:), nodes(j,:)));
            M(i,j)=phi(t,c);
    end
end
M(1:Nz, Nz+1) = nodes(:,1);
M(1:Nz, Nz+2) = nodes(:,2);
M(1:Nz, Nz+3) = nodes(:,3);
M(Nz+1, 1:Nz) = nodes(:,1)';
M(Nz+2, 1:Nz) = nodes(:,2)';
M(Nz+3, 1:Nz) = nodes(:,3)';
M(Nz+1:Nz+3,Nz+1:Nz+3)=zeros(3,3);

for i=1:Nz
    f(i) = func(nodes(i,:));
end
f(Nz+1:Nz+3)=zeros(3,1);

a = M\f;
end
%==================================================================
function vec = find_neighbours(node,u, nr)
  dist =real(acos(node*(u')));
  [sorted_dist,I] = sort(dist);
  vec = zeros(nr,3);
  vec(1:nr,:) = u(I(2:nr+1),:);
end
%==================================================================
function g = geodesic_dist(x,y)
g = acos(x*(y'));
end
%==================================================================
%Franke functions
function f = func(u)
x = u(1); y = u(2); z=u(3);
f = 0.75*exp(-((9*x-2)^2)/4 - ((9*y-2)^2)/4 - ((9*z-2)^2)/4)  ...
+ 0.75*exp(-((9*x+1)^2)/49 - ((9*y+1)^2)/10 - ((9*z+1)^2)/10) ...
+ 0.5*exp(-((9*x-7)^2)/4 - ((9*y-3)^2)/4 - ((9*z-5)^2)/4) ...
- 0.2*exp (-((9*x-4))^2 - (9*y-7)^2 - (9*z-5)^2); %f1
% f = (1.25+cos(5.4*y))*cos(6*z)/(6+6*(3*x-1)^2); %f2
% f = exp(-(81/16)*((x-0.5)^2+(y-0.5)^2+(z-0.5)^2))/3; %f3 
% f =  exp(-(81/4)*((x-0.5)^2+(y-0.5)^2+(z-0.5)^2))/3; %f4
end
%==================================================================
function f = phi(t,c)
g = 2 * sin(t/2);
f = 1 / sqrt(g^2 + c^2); %%inverse mq
end



