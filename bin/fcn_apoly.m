function [  ] = fcn_apoly( hObject, eventdata ,fig )
data=guidata(fig);
%idx=get(data.handles.listbox(1),'Value');

lenStep=length(data.userdata.step);
StepStr=data.userdata.step{lenStep};

D=StepStr.D;
T=StepStr.T;
X=StepStr.X;


lenStep=length(data.userdata.step);
StepStr=data.userdata.step{lenStep};

bs=str2num(get(data.handles.aPoly_edit(1),'String'));
eps_std=str2num(get(data.handles.aPoly_edit(2),'String'));
n=str2num(get(data.handles.aPoly_edit(3),'String'));
expand=str2num(get(data.handles.aPoly_edit(4),'String'));
peakW=str2num(get(data.handles.aPoly_edit(5),'String'));

%D=StepStr.D;
%T=StepStr.T;
X=StepStr.X;



opt.method='sim';
opt.maxST=eps_std;
opt.ExpLevel=expand;
opt.plot=0;

[R,opt]=calcReg(D,T,X,opt);

R(:,1:bs)=1;
R(:,end-bs+1:end)=1;

for i=1:size(D,1)
    if mod(i,10)==0
        disp([num2str(round(i/size(D,1)*100)) ' %'])
    end
    Xt=X(R(i,:));
    Dt=D(i,R(i,:));
    p = polyfit(Xt(:),Dt(:),n);
    z = polyval(p,X);
    D(i,:) = D(i,:)-z(:)';
end

StepStr.D=D;


[U,S,V]=svd(StepStr.D);
StepStr.U=U(:,1:data.userdata.maxSVD);
StepStr.S=S(1:data.userdata.maxSVD,1:data.userdata.maxSVD);
StepStr.V=V(:,1:data.userdata.maxSVD);


data.userdata.targets=[];
data.userdata.h_target=[];
data.userdata.seltar=0;

%Fill new one
data.userdata.step{lenStep+1}=StepStr;
Cstr=get(data.handles.listbox(1),'String');
Cstr{lenStep+1}=['aPoly(' get(data.handles.aPoly_edit(1),'String') ','...
    get(data.handles.aPoly_edit(2),'String') ','...
    get(data.handles.aPoly_edit(3),'String') ','...
    get(data.handles.aPoly_edit(4),'String') ','...
    get(data.handles.aPoly_edit(5),'String') ')'];

set(data.handles.listbox(1),'String',Cstr,'Value',lenStep+1);

guidata(fig,data);

updateGUI_PP(fig);


end





function [ R,opt ] = calcReg( D,T,X,opt )
T=T(:);
X=X(:)';
R=zeros(size(D));

switch opt.method
    case 'sim'
        Dx=zeros(size(D));
        %DDx=zeros(size(D));
        
        for i=1:length(T)
            Dx(i,1:end-1)=diff(D(i,:))./diff(X) ;
            %DDx(i,1:end-1)=diff(Dx(i,:))./diff(X) ;
        end
        Dx(:,end)=Dx(:,end-1);
        %DDx(:,end)=DDx(:,end-1);
        
        
        aDx=abs(Dx);
        Ex=aDx./max(max(aDx));
        
        maxST=opt.maxST;
        
        eps=0.3;
        b=1;
        while b
            RVal=Ex(Ex<=eps);
            me=mean(RVal);
            st=std(RVal);
            if sum(RVal<me-maxST*st)+sum(RVal>me+maxST*st)>0
                eps=eps/2;
            else
                b=0;
            end
        end
        
        ES=(Ex<=eps);
        PR=~ES;
        %Expand
        
        ExpLevel=opt.ExpLevel;
        ePR=PR;
        
        for i=1:size(PR,1)
            ePR(i,:)=smooth(ePR(i,:),2*ExpLevel+1);
            
        end
        for j=1:size(ES,2)
            ePR(:,j)=smooth(ePR(:,j),2*ExpLevel+1);
        end
        PR=(ePR>0);
        
        R=~PR;
        
        if opt.plot ==2
            figure();
            
            mesh(X,T,D)
            
            hold on
            [XX,YY]=meshgrid(X,T);
            plot3(XX(R),YY(R),D(R),'rx');
            
        elseif opt.plot ==1
            plot(X,D)
            hold on
            plot(X,Dx,'r')
        end
        
        
end
end

