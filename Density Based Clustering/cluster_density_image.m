function [C,cl] = cluster_density_image(J)
[m, n, p] = size(J);
data = reshape(double(J), m*n, p);
dist = pdist2(data,data);


percent = 2.0;
N = size(dist,1);
position = round((N*(N-1)/2)*percent/100);

tri_dist = triu(dist,1);
s_dist = sort(tri_dist(tri_dist~=0));
dc = s_dist(position);

fprintf('Computing Rho with gaussian kernel of radius: %12.6f\n', dc);


rho = sum(exp(-(dist./dc).^2),2);

max_dist=max(max(dist));

[rho_sorted,ordrho]=sort(rho,'descend');


delta(ordrho(1))=-1;

nneigh(ordrho(1))=0;

for k=2:N
   delta(ordrho(k))=max_dist;
   for p=1:k-1
     if(dist(ordrho(k),ordrho(p))<delta(ordrho(k)))
        delta(ordrho(k))=dist(ordrho(k),ordrho(p));
        nneigh(ordrho(k))=ordrho(p);   
     end
   end
end

delta(ordrho(1))=max(delta(:));


disp('Generated file:DECISION GRAPH')
figure
subplot(1,2,1)
plot(rho(:),delta(:),'o','MarkerSize',5,'MarkerFaceColor','b','MarkerEdgeColor','r');
title ('Decision Graph','FontSize',15.0)
xlabel ('\rho')
ylabel ('\delta')


for i=1:N
  gamma(i)=rho(i)*delta(i);
end
[gamma_sorted,ordgamma]=sort(gamma,'descend');
subplot(1,2,2)
plot(gamma_sorted(:),'o','MarkerSize',5,'MarkerFaceColor','b','MarkerEdgeColor','r');
title ('Gamma Graph','FontSize',15.0)


q = input('q=');


NCLUST=0;

for b=1:N
  cl(b)=-1;
end 

C = zeros(q-1,3);
for n=1:N
  if gamma(n) > gamma_sorted(q)
     NCLUST=NCLUST+1;
     cl(n)=NCLUST;
     C(NCLUST,:) = data(n,:);
     icl(NCLUST)=n;
  end
end
fprintf('NUMBER OF CLUSTERS: %i \n', NCLUST);


disp('assignation');
for i=1:N
  if (cl(ordrho(i))==-1)
    cl(ordrho(i))=cl(nneigh(ordrho(i)));
  end
end


for i=1:N
  halo(i)=cl(i);
end 
if (NCLUST>1)

  for i=1:NCLUST
    border_rho(i)=0;
  end

  for i=1:N-1
    for j=i+1:N

      if ((cl(i)~=cl(j))&& (dist(i,j)<=dc))
        rho_aver=(rho(i)+rho(j))/2.; 
        if (rho_aver>border_rho(cl(i)))
          border_rho(cl(i))=rho_aver;
        end
        if (rho_aver>border_rho(cl(j)))
          border_rho(cl(j))=rho_aver;
        end
      end
    end
  end

  for i=1:N
    if (rho(i)<border_rho(cl(i)))
      halo(i)=0;
    end
  end 
end


for i=1:NCLUST
  nc=0; 
  nh=0; 
  for j=1:N
    if (cl(j)==i)
      nc=nc+1;
    end
    if (halo(j)==i)
      nh=nh+1;   
    end
  end
  fprintf('CLUSTER: %i CENTER: %i ELEMENTS: %i CORE: %i HALO: %i \n', i,icl(i),nc,nh,nc-nh); 
end