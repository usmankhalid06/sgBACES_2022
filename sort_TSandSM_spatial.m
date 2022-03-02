function [srtd_Zt,srtd_Zs,ind_Zs]=sort_TSandSM_spatial(TC,SM,Zt,Zs,Kc,M,Ks,nC,R) 
   srcs = size(TC,2);
   V = size(Zs,2);
   for j=1:srcs
       if R ==1
           [~, ind_Zs(j)]  = max(abs(corr(abs(SM(j,:)'),abs(Zs'))));
           srtd_Zs(j,:) =  Zs(ind_Zs(j),:);
           srtd_Zt(:,j) = sign(corr(TC(:,j),Zt(:,ind_Zs(j))))*Zt(:,ind_Zs(j));
       else
           if j<=nC
               [~, ind_Zs(j)]  = max(abs(corr(abs(SM(j,:)'),abs(Zs(1:Kc,:)'))));
               srtd_Zs(j,:) =  Zs(ind_Zs(j),:);
               srtd_Zt(:,j) = sign(corr(TC(:,j),Zt(:,ind_Zs(j))))*Zt(:,ind_Zs(j));
           else
               [~, ind_Zs(j)]  = max(abs(corr(abs(SM(j,:)'),abs([zeros(Kc,V); Zs(Kc+1:Kc+M*Ks,:)]'))));
               srtd_Zs(j,:) =  Zs(ind_Zs(j),:);
               srtd_Zt(:,j) = sign(corr(TC(:,j),Zt(:,ind_Zs(j))))*Zt(:,ind_Zs(j));
           end
       end
   end  
end
