function [tao, FPI] = FMM_2(x, Perm)
    lx=80;ly=80;         % Dimension of the blocks
    nx=60;ny=220;        % number of blocks in x- and y-direction 
    phi=0.2;

    mio=0.88; %cp
    ct=1e-5; %  1/psi
    
    ns=size(x,2)/2;      % Number of Source points
    dx=lx/nx;
    dy=ly/ny;
    ngset=zeros(ny,nx);  % Set of Neighbour blocks
    fset=ones(ny,nx);    % Set of Frozen Blocks
    taoset=inf(ny,nx);     % Value of tao for blocks
    taocal=inf(ny,nx);
    tao=inf(ny,nx);  %  Matrix assisting finding value of tao and index
    taoind=inf(ny*nx,3);  %  Matrix assisting finding value of tao and index
    FPI=zeros(ny,nx);   % Frozen Points Index
    FPInd = zeros(ny,nx);
    for jj=1:ny*nx;
       taoind(jj,1)=jj; 
    end
    minmat=zeros(ny*nx);
    arrangment=zeros(ny*nx,1);

    nwbx = x(1:2:end);
    nwby = x(2:2:end);

    alfa=sqrt(Perm./(3792.*mio.*ct.*phi));

    for i=1:ns
        fset((nwbx(i)-1)*ny+nwby(i))=0;
        taoset((nwbx(i)-1)*ny+nwby(i))=0;
        taocal((nwbx(i)-1)*ny+nwby(i))=0;
        tao((nwbx(i)-1)*ny+nwby(i))=0;
        FPI((nwbx(i)-1)*ny+nwby(i)) = i;
        FPInd((nwbx(i)-1)*ny+nwby(i)) = i;
    end

    fcs=2:ny-1;
    lcs=(nx-1)*ny+2:(ny)*nx-1;
    frs=(ny+1):ny:(nx-2)*ny+1;
    lrs=2*ny:ny:(nx-1)*ny;
    v=0;

    for rr=1:ns
        taoind(rr,2)=((nwbx(rr)-1)*ny+nwby(rr));
        taoind(rr,3)=0;
        arrangment(rr)=(nwbx(rr)-1)*ny+nwby(rr);
    end

  while sum(sum(fset))~=0 && v<ny*nx
        v=v+1;

      if v==1
        nh=ns;
      else
         nh=1;
      end

       aa=zeros(nh,1);

        for j=1:nh
            if v==1
                ctr=(v-1)*ns;
            else
                ctr=ns+(v-2);
            end
        end

    ngset;
    fset;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Determination of Tao
    nbm=zeros(4*nh,1);

                             %%% Specification of neighbours for tao
    for j=1:nh                         %%%%%%%%%%% calculation

        if taoind(ctr+j,2)==1

            if fset(taoind(ctr+j,2)+1)==1
                nbm(1+(j-1)*nh*2)=taoind(ctr+j,2)+1;
            end                                      %%%%%% First block

            if  fset(taoind(ctr+j,2)+nx)==1
                nbm(2+(j-1)*nh*2)=taoind(ctr+j,2)+ny;
            end


        elseif  taoind(ctr+j,2)==ny*nx

           if fset(taoind(ctr+j,2)-1)==1
             nbm(1+(j-1)*nh*2)=taoind(ctr+j,2)-1;
           end                                      %%%%%% ny*nx block

            if fset(taoind(ctr+j,2)-nx)==1
                nbm(2+(j-1)*nh*2)=taoind(ctr+j,2)-ny;
            end

        elseif taoind(ctr+j,2)==ny

           if fset(taoind(ctr+j,2)-1)==1
                nbm(1+(j-1)*nh*2)=taoind(ctr+j,2)-1;
           end                                                  %%%%%% ny block

            if fset(taoind(ctr+j,2)+ny)==1
                nbm(2+(j-1)*nh*2)=taoind(ctr+j,2)+ny;
            end

        elseif taoind(ctr+j,2)==(nx-1)*ny+1

            if fset(taoind(ctr+j,2)+1)==1
                nbm(1+(j-1)*nh*2)=taoind(ctr+j,2)+1;
            end                                      %%%%%% nx block

            if fset(taoind(ctr+j,2)-ny)==1
                  nbm(2+(j-1)*nh*2)=taoind(ctr+j,2)-ny;
            end


        elseif  size(find(fcs==taoind(ctr+j,2)),2)~=0                      %%%%%%% For first coulmn blocks

              if fset((taoind(ctr+j,2))+1)==1
                nbm(1+(j-1)*nh*2)=taoind(ctr+j,2)+1;
              end

              if fset((taoind(ctr+j,2))-1)==1  
                nbm(2+(j-1)*ns*2)=taoind(ctr+j,2)-1;
              end


              if fset((taoind(ctr+j,2))+ny)==1
                nbm(4+(j-1)*nh*2)=taoind(ctr+j,2)+ny;
              end


         elseif  size(find(lcs==taoind(ctr+j,2)),2)~=0                      %%%%%%% For last coulmn blocks

              if fset((taoind(ctr+j,2))+1)==1
                nbm(1+(j-1)*nh*2)=taoind(ctr+j,2)+1;
              end

                if fset((taoind(ctr+j,2))-1)==1  
                    nbm(2+(j-1)*nh*2)=taoind(ctr+j,2)-1;
                end


                if fset((taoind(ctr+j,2))-ny)==1
                    nbm(4+(j-1)*ns*2)=taoind(ctr+j,2)-ny;
                end   


             elseif  size(find(frs==taoind(ctr+j,2)),2)~=0                      %%%%%%% For first row set blocks

                if fset((taoind(ctr+j,2))+1)==1
                    nbm(1+(j-1)*nh*2)=taoind(ctr+j,2)+1;
                end

                if fset((taoind(ctr+j,2))+ny)==1   
                    nbm(2+(j-1)*nh*2)=taoind(ctr+j,2)+ny;
                 end


                if fset((taoind(ctr+j,2))-ny)==1
                    nbm(4+(j-1)*nh*2)=taoind(ctr+j,2)-ny;
                end   


            elseif  size(find(lrs==taoind(ctr+j,2)),2)~=0                      %%%%%%% For last row set blocks

               if fset((taoind(ctr+j,2))-1)==1
                nbm(1+(j-1)*nh*2)=taoind(ctr+j,2)-1;
               end

              if fset((taoind(ctr+j,2))+ny)==1   
                nbm(2+(j-1)*nh*2)=taoind(ctr+j,2)+ny;
              end


            if fset((taoind(ctr+j,2))-ny)==1
                nbm(4+(j-1)*nh*2)=taoind(ctr+j,2)-ny;
            end   



        else        %%%%%%%%%%%%%%% neighbour specification for middle blocks

            if fset((taoind(ctr+j,2))+1)==1
                nbm(1+(j-1)*nh*2)=taoind(ctr+j,2)+1;
            end

            if fset((taoind(ctr+j,2))-1)==1  %for middle blocks
                nbm(2+(j-1)*nh*2)=taoind(ctr+j,2)-1;
            end

            if fset((taoind(ctr+j,2))-ny)==1
                nbm(3+(j-1)*nh*2)=taoind(ctr+j,2)-ny;
            end

            if fset((taoind(ctr+j,2))+ny)==1
                nbm(4+(j-1)*nh*2)=taoind(ctr+j,2)+ny;
            end


        end  
    end
        %%%%%%%%%%%%%%%%%%%%% End of neighbour specification

    nbm=unique(nbm);
    nnbm=size(nonzeros(nbm),1);
    nbm=nonzeros(nbm);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% START OF TAO CALCULATION
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for jj=1:nnbm
        if nbm(jj)==1                    % Calculation of tao for first block(if it is a neighbour)

        taoyp = taocal(nbm(jj) + 1);
        taoxp = taocal(nbm(jj) + ny);
        mat = [taoxp taoyp];
        [mat, n] = sort(mat);

        switch n(1)
            case 1
                FPI(nbm(jj)) = FPI(nbm(jj) + ny);
            case 2
                FPI(nbm(jj)) = FPI(nbm(jj) + 1);
            otherwise
                break
        end

        taox(nbm(jj))=taocal(1+ny);           
        taoy(nbm(jj))=taocal(2);

           if taox(nbm(jj))==inf
              tao(nbm(jj))=taoy(nbm(jj))+(dy/(alfa(nbm(jj))));
                 taoset(nbm(jj))=tao(nbm(jj));
         elseif taoy(nbm(jj))==inf
                   tao(nbm(jj))=taox(nbm(jj))+(dx/(alfa(nbm(jj))));
                      taoset(nbm(jj))=tao(nbm(jj));
          else
                    a=taox(nbm(jj));
                    b=dx;
                    c=taoy(nbm(jj));
                    d=dy;  
                    e=1/alfa(nbm(jj));
                    tao(nbm(jj))= (a*d^2 + b^2*c + b*d*(- a^2 + 2*a*c + b^2*e^2 - c^2 + d^2*e^2)^(1/2))/(b^2 + d^2);
                     taoset(nbm(jj))=tao(nbm(jj));
                     if isreal( tao(nbm(jj)))~=1
                     min_complex=min(a,c);
                         if min_complex==a
                         tao(nbm(jj))=taox(nbm(jj))+(dx/(alfa(nbm(jj))));
                         else
                         tao(nbm(jj))=taoy(nbm(jj))+(dy/(alfa(nbm(jj))));
                         end
                     end
           end

        elseif nbm(jj)==ny;              % Calculation of tao for last block in first coulmn(if it is a neighbour)
            taoym = taocal(nbm(jj) - 1);
            taoxp = taocal(nbm(jj) + ny);
            mat = [taoxp taoym];
            [mat, n] = sort(mat);

            switch n(1)
                case 1
                    FPI(nbm(jj)) = FPI(nbm(jj) + ny);
                case 2
                    FPI(nbm(jj)) = FPI(nbm(jj) - 1);
                otherwise
                    break
            end

          taox(nbm(jj))=taocal(2*ny);           
          taoy(nbm(jj))=taocal(ny-1);

                if taox(nbm(jj))==inf
                 tao(nbm(jj))=taoy(nbm(jj))+(dy/(alfa(nbm(jj))));
                 taoset(nbm(jj))=tao(nbm(jj));
                elseif taoy(nbm(jj))==inf
                   tao(nbm(jj))=taox(nbm(jj))+(dx/(alfa(nbm(jj))));
                   taoset(nbm(jj))=tao(nbm(jj));
                else
                    a=taox(nbm(jj));
                    b=dx;
                    c=taoy(nbm(jj));
                    d=dy;  
                    e=1/alfa(nbm(jj));
                    tao(nbm(jj))= (a*d^2 + b^2*c + b*d*(- a^2 + 2*a*c + b^2*e^2 - c^2 + d^2*e^2)^(1/2))/(b^2 + d^2);
                     taoset(nbm(jj))=tao(nbm(jj)); 
                if isreal( tao(nbm(jj)))~=1
                     min_complex=min(a,c);
                         if min_complex==a
                         tao(nbm(jj))=taox(nbm(jj))+(dx/(alfa(nbm(jj))));
                         else
                         tao(nbm(jj))=taoy(nbm(jj))+(dy/(alfa(nbm(jj))));
                         end
                end


                end



        elseif nbm(jj)==(nx-1)*ny+1;    % Calculation of tao for (ny-1)*nx+1 block (first block in last column)
            taoyp = taocal(nbm(jj) + 1);
            taoxm = taocal(nbm(jj) - ny);
            mat = [taoxm taoyp];
            [mat, n] = sort(mat);

            switch n(1)
                case 1
                    FPI(nbm(jj)) = FPI(nbm(jj) - ny);
                case 2
                    FPI(nbm(jj)) = FPI(nbm(jj) + 1);
                otherwise
                    break
            end

             taox(nbm(jj))=taocal((nx-2)*ny+1);           
             taoy(nbm(jj))=taocal((nx-1)*ny+2);
                if taox(nbm(jj))==inf
                 tao(nbm(jj))=taoy(nbm(jj))+(dy/(alfa(nbm(jj))));
                 taoset(nbm(jj))=tao(nbm(jj));
                elseif taoy(nbm(jj))==inf
                   tao(nbm(jj))=taox(nbm(jj))+(dx/(alfa(nbm(jj)))); 
                   taoset(nbm(jj))=tao(nbm(jj));
                else
                    a=taox(nbm(jj));
                    b=dx;
                    c=taoy(nbm(jj));
                    d=dy;  
                    e=1/alfa(nbm(jj));
                    tao(nbm(jj))= (a*d^2 + b^2*c + b*d*(- a^2 + 2*a*c + b^2*e^2 - c^2 + d^2*e^2)^(1/2))/(b^2 + d^2);
                     taoset(nbm(jj))=tao(nbm(jj));                                 

                   if isreal( tao(nbm(jj)))~=1
                     min_complex=min(a,c);
                         if min_complex==a
                         tao(nbm(jj))=taox(nbm(jj))+(dx/(alfa(nbm(jj))));
                         else
                         tao(nbm(jj))=taoy(nbm(jj))+(dy/(alfa(nbm(jj))));
                         end
                   end
                end

        elseif  nbm(jj)==(ny)*(nx);   % For (ny)*(nx) block (last block in last column)
            taoym = taocal(nbm(jj) - 1);
            taoxm= taocal(nbm(jj) - ny);
            mat = [taoxm taoym];
            [mat, n] = sort(mat);

            switch n(1)
                case 1
                    FPI(nbm(jj)) = FPI(nbm(jj) - ny);
                case 2
                    FPI(nbm(jj)) = FPI(nbm(jj) - 1);
                otherwise
                    break
            end
            taox(nbm(jj))=taocal((nx-1)*ny);           
            taoy(nbm(jj))=taocal(nx*ny-1);
            if( taox(nbm(jj))==inf)
                tao(nbm(jj))=taoy(nbm(jj))+(dy/(alfa(nbm(jj))));
                taoset(nbm(jj))=tao(nbm(jj));
            elseif taoy(nbm(jj))==inf
                tao(nbm(jj))=taox(nbm(jj))+(dx/(alfa(nbm(jj)))); 
                taoset(nbm(jj))=tao(nbm(jj));
            else
                a=taox(nbm(jj));
                b=dx;
                c=taoy(nbm(jj));
                d=dy;  
                e=1/alfa(nbm(jj));
                tao(nbm(jj))= (a*d^2 + b^2*c + b*d*(- a^2 + 2*a*c + b^2*e^2 - c^2 + d^2*e^2)^(1/2))/(b^2 + d^2);
                taoset(nbm(jj))=tao(nbm(jj)); 

                if isreal( tao(nbm(jj)))~=1
                    min_complex=min(a,c);
                    if min_complex==a
                        tao(nbm(jj))=taox(nbm(jj))+(dx/(alfa(nbm(jj))));
                    else
                         tao(nbm(jj))=taoy(nbm(jj))+(dy/(alfa(nbm(jj))));
                    end
                end
            end


        elseif size(find(fcs==nbm(jj)),2)~=0         % Tao calculation for first column set (fcs)
        taoyp = taocal(nbm(jj) + 1);
        taoym = taocal(nbm(jj) - 1);
        taoxp = taocal(nbm(jj) + ny);
        mat = [taoxp taoyp taoym];
        [mat, n] = sort(mat);

        switch n(1)
            case 1
                FPI(nbm(jj)) = FPI(nbm(jj) + ny);
            case 2
                FPI(nbm(jj)) = FPI(nbm(jj) + 1);
            case 3
                FPI(nbm(jj)) = FPI(nbm(jj) - 1);
            otherwise
                break
        end

        taoy(nbm(jj))=min(taoym,taoyp);           
        taox(nbm(jj))=taocal(nbm(jj)+ny);
                if taox(nbm(jj))==inf
                 tao(nbm(jj))=taoy(nbm(jj))+(dy/(alfa(nbm(jj))));
                 taoset(nbm(jj))=tao(nbm(jj));
                elseif taoy(nbm(jj))==inf
                   tao(nbm(jj))=taox(nbm(jj))+(dx/(alfa(nbm(jj)))); 
                   taoset(nbm(jj))=tao(nbm(jj));
                else
                    a=taox(nbm(jj));
                    b=dx;
                    c=taoy(nbm(jj));
                    d=dy;  
                    e=1/alfa(nbm(jj));
                    tao(nbm(jj))= (a*d^2 + b^2*c + b*d*(- a^2 + 2*a*c + b^2*e^2 - c^2 + d^2*e^2)^(1/2))/(b^2 + d^2);
                    taoset(nbm(jj))=tao(nbm(jj)); 

                     if isreal( tao(nbm(jj)))~=1
                        min_complex=min(a,c);
                         if min_complex==a
                         tao(nbm(jj))=taox(nbm(jj))+(dx/(alfa(nbm(jj))));
                         else
                         tao(nbm(jj))=taoy(nbm(jj))+(dy/(alfa(nbm(jj))));
                         end
                     end

                end                                           % end of first column set

        elseif size(find(lcs==nbm(jj)),2)~=0         % Tao calculation for last column set (lcs)
        taoyp = taocal(nbm(jj) + 1);
        taoym = taocal(nbm(jj) - 1);
        taoxm = taocal(nbm(jj) - ny);
        mat = [taoxm taoyp taoym];
        [mat, n] = sort(mat);

        switch n(1)
            case 1
                FPI(nbm(jj)) = FPI(nbm(jj) - ny);
            case 2
                FPI(nbm(jj)) = FPI(nbm(jj) + 1);
            case 3
                FPI(nbm(jj)) = FPI(nbm(jj) - 1);
            otherwise
                break
        end

        taoy(nbm(jj))=min(taoym,taoyp);           
        taox(nbm(jj))=taocal(nbm(jj)-ny);
                if taox(nbm(jj))==inf
                 tao(nbm(jj))=taoy(nbm(jj))+(dy/(alfa(nbm(jj))));
                 taoset(nbm(jj))=tao(nbm(jj));
                elseif taoy(nbm(jj))==inf
                   tao(nbm(jj))=taox(nbm(jj))+(dx/(alfa(nbm(jj))));
                   taoset(nbm(jj))=tao(nbm(jj));
                else
                    a=taox(nbm(jj));
                    b=dx;
                    c=taoy(nbm(jj));
                    d=dy;  
                    e=1/alfa(nbm(jj));
                    tao(nbm(jj))= (a*d^2 + b^2*c + b*d*(- a^2 + 2*a*c + b^2*e^2 - c^2 + d^2*e^2)^(1/2))/(b^2 + d^2);
                    taoset(nbm(jj))=tao(nbm(jj)); 


                      if isreal( tao(nbm(jj)))~=1
                          min_complex=min(a,c);
                         if min_complex==a
                         tao(nbm(jj))=taox(nbm(jj))+(dx/(alfa(nbm(jj))));
                         else
                         tao(nbm(jj))=taoy(nbm(jj))+(dy/(alfa(nbm(jj))));
                         end
                     end

                 end                                           
                                                  % end of last column set


        elseif size(find(frs==nbm(jj)),2)~=0         % Tao calculation for first row set (frs)
        taoyp = taocal(nbm(jj) + 1);
        taoxp = taocal(nbm(jj) + ny);
        taoxm = taocal(nbm(jj) - ny);
        mat = [taoxp taoxm taoyp];
        [mat, n] = sort(mat);

        switch n(1)
            case 1
                FPI(nbm(jj)) = FPI(nbm(jj) + ny);
            case 2
                FPI(nbm(jj)) = FPI(nbm(jj) - ny);
            case 3
                FPI(nbm(jj)) = FPI(nbm(jj) + 1);
            otherwise
                break
        end

        taoy(nbm(jj))=taoyp;           
        taox(nbm(jj))=min(taoxp, taoxm);
                if taox(nbm(jj))==inf
                 tao(nbm(jj))=taoy(nbm(jj))+(dy/(alfa(nbm(jj))));
                 taoset(nbm(jj))=tao(nbm(jj));
                elseif taoy(nbm(jj))==inf
                   tao(nbm(jj))=taox(nbm(jj))+(dx/(alfa(nbm(jj))));
                   taoset(nbm(jj))=tao(nbm(jj));
                else
                    a=taox(nbm(jj));
                    b=dx;
                    c=taoy(nbm(jj));
                    d=dy;  
                    e=1/alfa(nbm(jj));
                    tao(nbm(jj))= (a*d^2 + b^2*c + b*d*(- a^2 + 2*a*c + b^2*e^2 - c^2 + d^2*e^2)^(1/2))/(b^2 + d^2);
                    taoset(nbm(jj))=tao(nbm(jj));   

                     if isreal( tao(nbm(jj)))~=1
                         min_complex=min(a,c);
                         if min_complex==a
                         tao(nbm(jj))=taox(nbm(jj))+(dx/(alfa(nbm(jj))));
                         else
                         tao(nbm(jj))=taoy(nbm(jj))+(dy/(alfa(nbm(jj))));
                         end
                     end

                 end                                           
                                                   % end of first row set  



         elseif size(find(lrs==nbm(jj)),2)~=0         % Tao calculation for last row set (frs)
            taoym = taocal(nbm(jj) -1);
            taoxp = taocal(nbm(jj) + ny);
            taoxm = taocal(nbm(jj) - ny);
            mat = [taoxp taoxm taoym];
            [mat, n] = sort(mat);

            switch n(1)
                case 1
                    FPI(nbm(jj)) = FPI(nbm(jj) + ny);
                case 2
                    FPI(nbm(jj)) = FPI(nbm(jj) - ny);
                case 3
                    FPI(nbm(jj)) = FPI(nbm(jj) - 1);
                otherwise
                    break
            end

            taoy(nbm(jj))=taoym;           
            taox(nbm(jj))=min(taoxp, taoxm);
                if taox(nbm(jj))==inf
                 tao(nbm(jj))=taoy(nbm(jj))+(dy/(alfa(nbm(jj))));
                 taoset(nbm(jj))=tao(nbm(jj));
            elseif taoy(nbm(jj))==inf
                   tao(nbm(jj))=taox(nbm(jj))+(dx/(alfa(nbm(jj))));
                   taoset(nbm(jj))=tao(nbm(jj));
             else
                    a=taox(nbm(jj));
                    b=dx;
                    c=taoy(nbm(jj));
                    d=dy;  
                    e=1/alfa(nbm(jj));
                    tao(nbm(jj))= (a*d^2 + b^2*c + b*d*(- a^2 + 2*a*c + b^2*e^2 - c^2 + d^2*e^2)^(1/2))/(b^2 + d^2);
                    taoset(nbm(jj))=tao(nbm(jj));    


                     if isreal( tao(nbm(jj)))~=1
                          min_complex=min(a,c);
                         if min_complex==a
                         tao(nbm(jj))=taox(nbm(jj))+(dx/(alfa(nbm(jj))));
                         else
                         tao(nbm(jj))=taoy(nbm(jj))+(dy/(alfa(nbm(jj))));
                         end
                     end 
               end       
                                                   % end of last row set  
        else                                     % For inside (middle) blocks
            taoym = taocal(nbm(jj) - 1);
            taoxp = taocal(nbm(jj) + ny);
            taoxm = taocal(nbm(jj) - ny);
            taoyp = taocal(nbm(jj) + 1); 
            mat = [taoxp taoxm taoyp taoym];
            [mat, n] = sort(mat);

            switch n(1)
                case 1
                    FPI(nbm(jj)) = FPI(nbm(jj) + ny);
                case 2
                    FPI(nbm(jj)) = FPI(nbm(jj) - ny);
                case 3
                    FPI(nbm(jj)) = FPI(nbm(jj) + 1);
                case 4
                    FPI(nbm(jj)) = FPI(nbm(jj) - 1);
                otherwise
                    break
            end
         taox(nbm(jj))=min(taocal(nbm(jj)+ny),taocal(nbm(jj)-ny));        
         taoy(nbm(jj))=min(taocal(nbm(jj)-1),taocal(nbm(jj)+1));
                if taox(nbm(jj))==inf
                 tao(nbm(jj))=taoy(nbm(jj))+(dy/(alfa(nbm(jj))));
                 taoset(nbm(jj))=tao(nbm(jj));
                elseif taoy(nbm(jj))==inf
                   tao(nbm(jj))=taox(nbm(jj))+(dx/(alfa(nbm(jj))));
                   taoset(nbm(jj))=tao(nbm(jj));
                else
                    a=taox(nbm(jj));
                    b=dx;
                    c=taoy(nbm(jj));
                    d=dy;  
                    e=1/alfa(nbm(jj));
                    tao(nbm(jj))= (a*d^2 + b^2*c + b*d*(- a^2 + 2*a*c + b^2*e^2 - c^2 + d^2*e^2)^(1/2))/(b^2 + d^2);
                    taoset(nbm(jj))=tao(nbm(jj)); 

                    if isreal( tao(nbm(jj)))~=1
                          min_complex=min(a,c);
                         if min_complex==a
                         tao(nbm(jj))=taox(nbm(jj))+(dx/(alfa(nbm(jj))));
                         else
                         tao(nbm(jj))=taoy(nbm(jj))+(dy/(alfa(nbm(jj))));
                         end
                     end

                end  
                                                 % End of middle block
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% End of Tao
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Calculation


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Finding the
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% minimum value
    minmat=taoset;
    minmat(taoset==0)=nan;
    [M,I]=min(minmat(:));
    if isnan(M)
        break
    end
    taoset;
    fset(I)=0;              %%%%%%%%%%%%%% Setting the block to frozen
    Col = (I - mod(I,ny))/220 + 1;
    Row = mod(I,ny);
    if(Row==0)
        Row=220;
    end
    FPInd(I) = FPI(Row, Col);
    ngset(I)=0;          %%%%%%%%%%%%%% 
    taoset(I)=0 ;         %%%%%%%%%%%%%%  Setting tao to zero
    taocal(I)=M;
    tao(I)=M;
    arrangment(ns+v)=I;
    taoind(ns+v,2)=I;
    taoind(ns+v,3)=M;
    taoind;

end

