function [DelpT,DelpReal, TFL]=FMM(x)
    nx=25;ny=25;nz=1;

    global NFE
    NFE=NFE+1;
 
    x(1:2:end)=min(floor(1+(nx).*x(1:2:end)),(nx));
    x(2:2:end)=min(floor(1+(ny).*x(2:2:end)),(ny));
%     x = [15 15 6 20 10 4];

    %.... Number of wells
    nwd=numel(x)/2;

    %...Import Source Points 3D
%     for i=1:nwd
%         for j = 1:nz
%             WellLocation(nz*(i-1)+j,:)=[x(2*i-1:2*i) j];
%         end
%     end

    % 2 dimensinal
    WellLocation=zeros(nwd,2);
    for i=1:nwd
        WellLocation(i,:)=x(2*i-1:2*i);
    end


    %... Reservoir Data
    ct=5E-6;                      %.... Total Compressibility Value
    mio=0.88;                     %.... Viscosity
    dx=50;dy=50;                  %.... Block Length
    Bo=1.12;                      %... FVF
    qwell=50;                     %... Production Flow Rate
    h=10;                         %... Reservoir Thickness
    rw=0.66;                      %... Well Radius
    % P_ResInitial=4407;          %.... Initial reservoir Pressure

    %... Import Perm and Porocity
    Perm=load('Perm-het-syn.txt');
    Perm=reshape(Perm,[ny nx]);
%     perm = zeros(nx,ny,nz);
%     for i=1:nz
%         perm(:,:,i) = Perm(:,:,i)';
%     end
%     Perm = perm;
%     clear perm
%     Perm=ones(nx,ny,nz)*20;
    % Poro=load('Poro.txt');
    % Poro=reshape(Poro,[nx ny]);
    Poro=0.15.*ones(ny,nx,nz);
    % Perm=5.*ones(nx,ny);

    %... Calculation of Hydraulic Diffusivity Value
    HDV=(1/dx)*sqrt(Perm./(3792.*mio*ct.*Poro));

    %.... Run Kroon Code to find time of flight
    TFL=msfm(HDV, WellLocation', true, false);

    %......
    [WellLocation,~,c]=unique(WellLocation,'rows','stable');
    %....

    nw=size(WellLocation,1);
    %... Drained Block Construction
    DBM=FDKK(nw,TFL);
    disp('Time to estimate the drainage boundary. (DBM): ');
    toc
    imagesc(DBM);
    title('Drainage boundary of production wells at corners')
    axis image
    %.... Showing Drained Blocks
    DBF=DBM;

    for i=1:nw
        DBF(WellLocation(i,1),WellLocation(i,2))=0;
    end

    %... Defining Well Index
    ro=0.14*sqrt(dx^2+dy^2);
    WI=(2*pi*0.001127*h)/mio/(log(ro/rw)).*Perm;

    %.... Blocks Pore volume
    BPV=dx*dy*h.*Poro;

    Well=PAM(nw,TFL,DBM,BPV);

    Delp=zeros(nw,1);
    DelpT=1;
    for i=1:nw

    %... Creating the matrix [TOF BlockPoreVolume]
    % W1=VPAM(i,TFL,DBM,BPV);
    VPAM=Well(i).No;
    int=Matrix(VPAM);

    %... Reservoir to Block Pressure Drop -- P_Reservoir - P_Block
    DelpResBlock=(qwell*Bo/ct)*(5.615/24).*(int)';

    %... Block to Well Pressure Drop -- P_Block - P_Well
    DelpBlockWell=qwell/WI(WellLocation(i,1),WellLocation(i,2));

    %... Reservoir to Well Pressure Drop -- P_Reservoir - P_Well
    Delp(i)=DelpResBlock+DelpBlockWell;
    % a=max(P_ResInitial-Delp(i),1);
    % BHP=(a)*BHP;
    DelpT=DelpT.*Delp(i);
    end

    DelpT=(DelpT)^(1/nw);
    DelpReal=zeros(nwd,1);
    for i=1:nwd
        if i==c(i)
    DelpReal(i)=Delp(c(i));
        end
    end
    %... BHP pressure
    % BHPi=P_ResInitial-DelpReal;
    % Delpi=DelpReal;
    % toc(Overal_time);
    % figure;
    % imagesc(TFL);
    % figure;
    % imagesc(DBF);
    % toc;
 end

function Well=PAM(nw,TFL,DBM,BPV)
    Value=zeros(nw,1);
    for i=1:nw
        Well(i).No=[sort(TFL(DBM==i)), BPV(DBM==i)];
        Value(i)=Well(i).No(end,1);
    end

    A=0.7.*min(Value);

    for i=1:nw
        d=size(Well(i).No,1);
        while eps>0
        if Well(i).No(d,1)<=A 
            break;
        end
        d=d-1;
        end
        Well(i).No(d+1:end,:)=[];
    end
end


function int=Matrix(VPAM)

    % VPAM(all(VPAM,2)==0,:)=[];
    VPAM(:,2)=cumsum(VPAM(:,2));

    %... eliminating Repetetive Elements
    RI1=all(diff(VPAM),2);
    VPAM(RI1==0,:)=[];

    %... Convertto real time
    VPAM(:,1)=((VPAM(:,1)).^2)/4;
    VPAM(:,2)=(1./VPAM(:,2));

    int=trapz(VPAM(:,1),VPAM(:,2));

end

