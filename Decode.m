function [findMem,Section]= Decode(Position,nX,bit,Dt,nXsum,Ex)   
%% ARRANGE
global ng_col ng_be
%Ex         =[ 1         2x     3 4 5  6  7  8  9  10 11 12 13 13 15 16 17
%             Index SectionName A d tw bf tf Ix Zx Sx rx Iy Zy Sy ry J Cw]
%Position=zeros(1,nX(1)*bit(1)+nX(2)*bit(2)+nX(3)*bit(3))
ncol=Dt.ncol; nbe=Dt.nbe; nMem=Dt.nMem; Nc=0; nbei=0; nbri=0; nstory=Dt.nstory; nbay=Dt.nbay; 
Ex(:,2)=[]; Ex1=[]; Ex2=[]; Ex3=[];
k1=size(Ex,1);
Ex1=Ex(1 : 2^bit(2),:); %beam
Ex2=Ex([169:234,268:329],:); %col
Ex3=Ex([235:267,268:298],:); %bracing
%EX1 contains properties of I-section, EX2 contains properties of H-section, 
%RealX contains All unknown variables of each particle in each iteration
%Possible divisoion (F1: [1 3], F2: [1 2 3 6], F3: [1 2 3 6 9])
findMem=[]; Section=[];
%TEST POSITION
     % nX=[3 3 9] bit=[7 7 5]
    % Position=[0 0 1 1 1 1  0 0 0 0 1 1 0 1 0 1 0 1 0 1 1 0 1 0 0 1 1 0 0 1 1 0 0 0 1 1 1 1 1 1 1  0 0 0 0 1 1 0 1 0 1 0 1 0 1 1 0 1 0 0 1 0 1 0 1 1 1 0 1 1 0 1 0 1 0 0 1 0 1 0 1 1 1 0 1 1 0 ];
%% CLASSIFICATION BEAMS & COLUMN INTO GROUPS
%For columns and beams with more than one group
ic=[1:ncol]; ic=reshape(ic,nstory,nbay+1);
GC=mat2cell(ic,ng_col,nbay+1);

ib=[ncol+1:ncol+nbe]; ib=reshape(ib,nstory,nbay);
GB=mat2cell(ib,ng_be,nbay);

%For bracing (divided X-bracing into 4 parts)
ibr=[ncol+nbe+1:ncol+nbe+nstory*nbay*4]';
GBr=mat2cell(ibr,4*ones(1,nstory*nbay),1);
%% ENCODING BINARY CODE TO SECTION'S INDEX
%nX: No. group of each types of member
%bit: No. bit per group
%position=nX*bit
%CLASSIFIED POSITION TO EACH TYPES
PC=Position(1:nX(1)*bit(1));
PBe=Position(length(PC)+1:length(PC)+nX(2)*bit(2));
PBr=Position(length(PC)+length(PBe)+1:end);
%CLASSIFIED POSITION FOR EACH TYPE TO GROUP
PC=mat2cell(PC,1,bit(1)*ones(1,nX(1)));
PBe=mat2cell(PBe,1,bit(2)*ones(1,nX(2)));
PBr=mat2cell(PBr,1,bit(3)*ones(1,nX(3)));
%index 0:bit-1
Nc=[0:bit(1)-1];
Nbe=[0:bit(2)-1];
Nbr=[0:bit(3)-2];
%% ENCODING BINARY CODE TO SECTION'S INDEX
%col
for i=1:nX(1)
    %Encoding
    Vc(i)=sum(PC{i}.*2.^Nc)+1; 
    %Define section properties
    gc=reshape(GC{i},size(GC{i},1)*size(GC{i},2),1);
    Section(gc,:) = repmat(Ex2(Vc(i),:),length(gc),1);
end
%beam
for i=1:nX(2)
    %Encoding
    Vbe(i)=sum(PBe{i}.*2.^Nbe)+1; 
    %Define section properties
    gb=reshape(GB{i},size(GB{i},1)*size(GB{i},2),1);
    Section(gb,:) = repmat(Ex1(Vbe(i),:),length(gb),1);
end
%brace
for i=1:nX(3)
    %Encoding
    gbr=reshape(GBr{i},size(GBr{i},1)*size(GBr{i},2),1);
    if PBr{i}(1)~=0 %First bit controls topology
        Vbr(i)=sum(PBr{i}(2:end).*2.^Nbr)+1; 
        %Define section properties
        Section(gbr,:) = repmat(Ex3(Vbr(i),:),length(gbr),1);
    else
        Section(gbr,:) = 0; %If no bracing, Area is 0
    end
end

%% RESULT
%!!Find Index of brace
    findMem=find(Section(:,2)~=0); %not include no brace
    Section=Section(findMem,:);
