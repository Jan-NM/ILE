function mOrte = multiframepoints(Orte,distance)
%%Script that detects if a point is found in several following pictures
%  the values are combined and the accuracy is improved
%  the old Orte matrix is updated with only one new improved value for
%  every multiframe point.
%
%created by: Manfred Kirchgessner < Manfred.Kirchgessner@web.de>, 
%            Frederik Grüll <Frederik.Gruell@kip.uni-heidelberg.de>%
%%

mOrte=zeros(size(Orte));
maximg = max(Orte(:,9));
cnt = zeros(maximg,1);
maxcnt = 0;
mubar = waitbar(0,'Sorting points in partitions...please wait');
for sl = 1:maximg
     Orte_frame = Orte( Orte(:,9) == sl, : ); 
     cnt = size(Orte_frame,1);
     if(cnt>maxcnt)
         maxcnt=cnt;
     end;
end

frames = zeros(maxcnt,size(Orte,2),maximg);
for sl=1:maximg
     Orte_frame = Orte( Orte(:,9) == sl, : ); 
     cnt(sl,1) = size(Orte_frame,1);
     frames(1:cnt(sl,1),:,sl)=Orte_frame;
end

newcnt = 1;


for sl=1:maximg
       
       for i = 1:cnt(sl)
          cluster =  frames(i,:,sl);
          if(cluster(1,1)>0)
                pt_on = 0;
                found =1;
                my = cluster(1,2);
                mx = cluster(1,3);


                while (found==1 && (sl+pt_on) <maximg)

                  Orte_nextpic = frames(:,:,sl+pt_on);            
                  Orte_nextpic = Orte_nextpic( abs(Orte_nextpic(:,2)-my) < distance  , : );

                  if(pt_on == 0)
                      Orte_nextpic = Orte_nextpic( abs( Orte_nextpic(:,2) - my ) > 0  , : );
                  end

                  if (isempty(Orte_nextpic) )
                      if (pt_on==0)
                          pt_on = pt_on + 1;
                      else
                          found = 0; 
                      end
                  else
                      Orte_nextpic = Orte_nextpic( abs(Orte_nextpic(:,3)-mx) < distance , : );

                      if (isempty(Orte_nextpic) )
                          found = 0;
                      else                    
                          cluster =[cluster;Orte_nextpic];
                          pt_on = pt_on + 1;
                          %remove clusterpoints in current picture in newOrte
                          for ki=1:size(Orte_nextpic,1)
                              %get index in newOrte of removable Points
                              now_frame = Orte_nextpic(ki,9);
                              kill = find( round(frames(:,2,now_frame)*10000) == round(Orte_nextpic(ki,2)*10000),1);
                              
                              if(~isempty(kill))
                                frames( kill , 1 ,now_frame) = 0;   
                              end
                          end

                      end
                  end
                end

                Q  = cluster(:,8);
                my = cluster(:,2);
                mx = cluster(:,3);
                dy = cluster(:,4);
                dx = cluster(:,5);
                sy = cluster(:,6);
                sx = cluster(:,7);

                Qges=sum( cluster(:,8) );
                Qmax=sum( cluster(:,1) );

                newmy = sum(my.*Q)/Qges;
                newmx = sum(mx.*Q)/Qges;

                newdy = sqrt(sum((Q.*dy).^2)/Qges^2 );
                newdx = sqrt(sum((Q.*dx).^2)/Qges^2 );

                newsy = sqrt(sum((sy.*sy + my.*my).*Q)/Qges-newmy*newmy);
                newsx = sqrt(sum((sx.*sx + mx.*mx).*Q)/Qges-newmx*newmx);
                imgno = sl;

                if(size(cluster,2)>9) %old versions only have 9 columns
                  Qgesold= mean(cluster(:,10));
                  mOrte(newcnt,1:10) =[Qmax,newmy,newmx,newdy,newdx,newsy,newsx,Qges,imgno,Qgesold];
                  newcnt=newcnt+1;
                else
                  mOrte(newcnt,1:9) =[Qmax,newmy,newmx,newdy,newdx,newsy,newsx,Qges,imgno];
                  newcnt=newcnt+1;
                end
           end

       end
       if (mod(sl,maximg/80)==0)
          waitbar(sl/maximg,mubar,['Searching for multiframe signals...please wait ' num2str(sl) '/' num2str(maximg) ] );
       end

end
mOrte=mOrte(1:newcnt-1,:);
 
close(mubar);