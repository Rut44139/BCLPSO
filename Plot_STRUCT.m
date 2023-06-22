function p1=Plot_STRUCT(axis_edge,arrow,bt,Name)
D=bt.D; G=bt.G; SP=bt.SUP; F=D.SeL; Prop=bt.P;
nNode=D.nNode;Node=D.Node; nMem=D.nMem; nbay=D.nbay; nstory=D.nstory;
ncol=D.ncol; nbe=D.nbe; 
N_fx=SP.N_fx;
CMxi=G.CMxi; CMyi=G.CMyi; CMxj=G.CMxj; CMyj=G.CMyj; NDOF_store=G.NDOF_store;
%-----------------------Separate type of dof------------------------------%
%-----------------Original coordinate of each  node-----------------------%
Cdx= Node(:,1);
Cdy= Node(:,2);
min_x=min(Cdx);     max_x=max(Cdx);
min_y=min(Cdy);     max_y=max(Cdy);
%----------------------New coordinate of each node------------------------%
for i=1:(nstory+1)*(nbay+1)
    NDOF(i,:)=NDOF_store{i};
end
for i=(nstory+1)*(nbay+1)+1:nNode
    NDOF(i,:)=[NDOF_store{i} 0];
end
%%----------------------Plot Original structure-----------------------%%
%------------------------Plot original structure--------------------------%

x1=[];y1=[];
for i=1:nMem
    X1(:,i)=[CMxi(i,1) CMxj(i,1)];
    Y1(:,i)=[CMyi(i,1) CMyj(i,1)];
end
%%----------------Plot Support reaction of structure------------------%%
%1. STRUCTURE
figure(4)
grid on
p1=plot(X1,Y1,'b','LineWidth',1.5);
title('STRUCTURE INDEX')
xlabel('Length(m)')
ylabel('Height(m)')
axis([min_x-axis_edge max_x+axis_edge min_y-axis_edge max_y+axis_edge]);
set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
hold on
%-----------------------------Plot data point-----------------------------%
%-------------------------------Plot Index--------------------------------%
for i=1:nNode
    x1=Cdx(i); y1=Cdy(i);
    x2=Cdx(i); y2=Cdy(i);
    Sm=plot(x1,y1,'r.','MarkerSize', 25); %Dot    
end
br=ncol+nbe+1:nMem; k=1; m=0;
for i=1:nMem
    if i >= ncol+nbe+1 %Position (Brace)
        x1=sum(X1(:,i))/2;
        y1=sum(Y1(:,i))/2-0.5;
    else %Position (Column & Beam)
        x1=sum(X1(:,i))/2;
        y1=sum(Y1(:,i))/2;
    end
   % if i <=ncol+nbe %Name (Column & Beam)
        name=Name(Prop(i)+1); 
        tx=text(x1,y1+0.1, name , 'FontSize',11,...
           'HorizontalAlignment','center',...
           'VerticalAlignment', 'bottom');
    
   % end

end

%---------------------------Plot external force---------------------------%
figure(5)
grid on
p1=plot(X1,Y1,'b','LineWidth',1.5);
title('STRUCTURE')
xlabel('Length(m)')
ylabel('Height(m)')
axis([min_x-axis_edge max_x+axis_edge min_y-axis_edge max_y+axis_edge]);
set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
hold on
p1=plot(X1,Y1,'b','LineWidth',1.5);
fdof=find(F(:,1)~=0); %find dof with applied force
for i=1:length(fdof)
    ndof=fdof(i);
    [rfN,cfN]=find(ndof==NDOF); 
    x=Cdx(rfN);
    y=Cdy(rfN);
    f=double(F(ndof,1));
    tex=annotation('textarrow');
    tex.Parent = gca; 
    tex.Color='r';
    if cfN==1              
       tex.X= [x-arrow x]; tex.Y= [y y]; 
       tex.String={'Fx= ',num2str(f)};
    elseif cfN==2  
       tex.X= [x x]; tex.Y= [y-arrow y];
       tex.String={'Fy= ',num2str(f)};
    else 
       tex.X= [x x-arrow-0.1]; tex.Y= [y+arrow-0.2 y+arrow-0.6]; 
       tex.String={'M= ',num2str(f)};
       line=annotation('line');
       line.Parent=gca;
       line.X=[x x+arrow-0.1];
       line.Y=[y+arrow-0.2 y+arrow-0.6];
       line.Color='r';
    end
end
%--------------------Plot support type (only fixed)-----------------------%
for i=1:length(N_fx)
    nN=N_fx(i);
    x=Cdx(nN);
    y=Cdy(nN);
    Sup=annotation('rectangle');
    Sup.Parent = gca; 
    Sup.Position=[x-0.5 y-0.2 1 0.4];
    Sup.Color='g';
    Sup.LineWidth=1.5;
end