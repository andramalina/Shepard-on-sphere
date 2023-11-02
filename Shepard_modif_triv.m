function z = Shepard_modif_triv(x,nodes,f,Nw,p)
N=length(nodes);
dist = zeros(1,N);
sorted_dist = dist;
salt = dist;
numer = 0;
denomin = 0;
i0=1;
for k = 1:N
    dk = acos(x*(nodes(k,:)'));
    if dk ~= 0
        %%% calculul lui Rwk
        xx = nodes(k,:);
        dist = acos(xx*(nodes'));
        [sorted_dist,I] = sort(dist);
        salt(Nw+2:N) = sorted_dist(Nw+2:N) - sorted_dist(Nw+1:N-1);
        saltmax = max(salt);
        saltmin = min(salt);
        saltmed = (saltmax + saltmin)/2;
        saltsemnif = saltmed; %min(3*saltmin,saltmed);
        for i = Nw+2:N
            if salt(i) >= saltsemnif
                i0 = i;
                break
            end
        end
        Wk = 0;
        Rwk = dist(I(i0));
% Rwk=abs(max(nodes)-min(nodes))/2*sqrt(9/N);
        if Rwk > dk
            Wk = ((Rwk - dk)/Rwk/dk)^p;

        end
        numer = numer + Wk.*f(k);
        denomin = denomin + Wk;

    end
end
z = numer./denomin;
end