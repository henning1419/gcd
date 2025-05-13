function [] = KeyPressFigure( src,event,fig )
data=guidata(fig);
if strcmp(event.Key,'delete')
    if data.userdata.seltar>0
        data.userdata.targets(data.userdata.seltar,:)=[];
        delete(data.userdata.h_target(data.userdata.seltar));
        data.userdata.h_target(data.userdata.seltar)=[];
        data.userdata.seltar=0;
        
        tmp_tar=data.userdata.targets;
        Cstr=cell(size(tmp_tar,1),1);
        for i=1:size(tmp_tar,1)
            Cstr{i}=[num2str(i) ': ' num2str(tmp_tar(i,1),6) ', ' num2str(tmp_tar(i,2),6)];
        end
        %%%
        set(data.handles.tarlisthandles,'String',Cstr(:),'Value',1);
        if data.userdata.seltar>0
            set(data.handles.tarlisthandles,'Value',data.userdata.seltar);
        end
    end
end


guidata(fig,data);


end






