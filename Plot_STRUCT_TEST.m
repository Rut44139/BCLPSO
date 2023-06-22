function p1=Plot_STRUCT_TEST(D,G)
nNode=D.nNode;Node=D.Node; nMem=D.nMem; 
ncol=D.ncol; nbe=D.nbe; nbr=D.nbr;
CMxi=G.CMxi; CMyi=G.CMyi; CMxj=G.CMxj; CMyj=G.CMyj; MDOF_store=G.MDOF_store; NDOF_store=G.NDOF_store;
%-----------------------Separate type of dof------------------------------%
%-----------------Original coordinate of each  node-----------------------%
Cdx= Node(:,1);
Cdy= Node(:,2);
min_x=min(Cdx);     max_x=max(Cdx);
min_y=min(Cdy);     max_y=max(Cdy);
%----------------------New coordinate of each node------------------------%

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
%axis([min_x-axis_edge max_x+axis_edge min_y-axis_edge max_y+axis_edge]);
set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
hold on
%-----------------------------Plot data point-----------------------------%
%-------------------------------Plot Index--------------------------------%
for i=1:nNode
    x1=Cdx(i); y1=Cdy(i);
    x2=Cdx(i); y2=Cdy(i);
    Sm=plot(x1,y1,'r.','MarkerSize', 25); %Dot    
    tx=text(x1,y1+0.1, num2str(i) , 'FontSize',11,...
           'HorizontalAlignment','center',...
           'VerticalAlignment', 'bottom');
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
    if i <=ncol+nbe %Name (Column & Beam)
        tx=text(x1,y1+0.1, num2str(i) , 'FontSize',11,...
           'HorizontalAlignment','center',...
           'VerticalAlignment', 'bottom');
    else %Name (Brace)
        m=m+1; 
        if m==2
            k=k+2;
            m=0;
        else
            tx=text(x1,y1+0.1, num2str(i) , 'FontSize',12,...
           'HorizontalAlignment','center',...
           'VerticalAlignment', 'bottom');
        end      
    end

end

