function [srtd_Zt,srtd_Zs,ind_Zt]=sort_TSandSM_temporal(TC,SM,Zt,Zs,Kc,M,Ks,nC,R) 
   srcs = size(TC,2);
   N = size(Zt,1);
   for j=1:srcs
       if R ==1
           [~, ind_Zt(j)]  = max(abs(corr(TC(:,j),Zt)));
           srtd_Zs(j,:) = Zs(ind_Zt(j),:);
           srtd_Zt(:,j) = sign(corr(TC(:,j),Zt(:,ind_Zt(j))))*Zt(:,ind_Zt(j));
       else
           if j<=nC
               [~, ind_Zt(j)]  = max(abs(corr(TC(:,j),Zt(:,1:Kc))));
               srtd_Zs(j,:) = Zs(ind_Zt(j),:);
               srtd_Zt(:,j) = sign(corr(TC(:,j),Zt(:,ind_Zt(j))))*Zt(:,ind_Zt(j));
           else
               [~, ind_Zt(j)]  = max(abs(corr(TC(:,j),[zeros(N,Kc) Zt(:,Kc+1:Kc+M*Ks)])));
               srtd_Zs(j,:) = Zs(ind_Zt(j),:);
               srtd_Zt(:,j) = sign(corr(TC(:,j),Zt(:,ind_Zt(j))))*Zt(:,ind_Zt(j));
           end
       end
   end  
end




