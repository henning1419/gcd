function buttonTargetUp( hObject, eventdata ,obj)

set(obj.mainfig,'WindowButtonUpFcn',[],'WindowButtonMotionFcn',[])
set([obj.mainaxes;allchild(obj.mainaxes)],'ButtonDownFcn',{@buttonTargetDown,obj});

% tmp_tar=data.userdata.targets;
% tmp_tar_handles=data.userdata.h_target;
int=[obj.target.start obj.target.end];
if int(1)~=int(2)
    int=[min(int) max(int)];
%     if isempty(tmp_tar)
        obj.activetarget=int;
        obj.activeplot=obj.tarplot;
        
        %%% targetboxes
%         set(data.handles.mPoly_listbox(1),'String',{['1: ' num2str(tmp_tar(1),6) ', ' num2str(tmp_tar(2),6)]},'Value',1);
        
%     else
%         [idx1,isIn1]=getnextInt(tmp_tar,int(1));
%         [idx2,isIn2]=getprevInt(tmp_tar,int(2));
%         if isIn1
%             idx1=idx1-1;
%             newleft=tmp_tar(idx1,1);
%         else
%             newleft=int(1);
%         end
%         if isIn2
%             idx2=idx2+1;
%             newright=tmp_tar(idx2,2);
%         else
%             newright=int(2);
%         end
%         tmp_tar(idx1:idx2,:)=[];
%         tmp_tar=[tmp_tar(1:idx1-1,:) ; [newleft newright] ; tmp_tar(idx1:end,:)];
%         
%        
%         delete(tmp_tar_handles(idx1:idx2));
%         tmp_tar_handles(idx1:idx2)=[];
%         tmp_tar_handles=tmp_tar_handles(:);
%         tmp_tar_handles=[tmp_tar_handles(1:idx1-1) ; data.tmp.tarplot ; tmp_tar_handles(idx1:end)];
%         yy=get(data.handles.axesD,'Ylim');
%         set(data.tmp.tarplot,'Xdata',[newleft newright newright newleft],'Ydata',[yy(1) yy(1) yy(2) yy(2)],'Facecolor',data.guiInfo.tarcol)
%         
%     end
    
% else
%     delete(obj.tarplot);
%     idx=0;
%     for i=1:size(tmp_tar,1)
%         if int(1)>=tmp_tar(i,1) && int(1)<=tmp_tar(i,2)
%             idx=i;
%             break
%         end
%     end
%     data.userdata.seltar=idx;
%     set(tmp_tar_handles(:),'Facecolor',data.guiInfo.tarcol)
%     if idx>0 
%         set(tmp_tar_handles(idx),'Facecolor',data.guiInfo.tarselcol)  
%     end
    
end




set([obj.mainaxes;allchild(obj.mainaxes)],'ButtonDownFcn',{@buttonTargetDown,obj});
% set(fig,'KeyPressFcn',{@KeyPressFigure,fig})
% disp(data.userdata.targets);
obj.target=rmfield(obj.target,'start');
obj.target=rmfield(obj.target,'end');
end



% 
% 
% 
% function [idx,isInInt]=getprevInt(tar,x)
% % idx = idx of prev Int \in [0 size(tar,1)]
% % isInInt = 0 - no, 1 - yes with index idx1+1
% if x<=tar(1,2)
%     idx=0;
%     if x>=tar(1,1)
%         isInInt=1;
%     else
%         isInInt=0;
%     end
%     return;
% end
% if x>tar(end,2)
%     idx=size(tar,1);
%     isInInt=0;
%     return;
% end
% for i=2:size(tar,1)
%     if x>tar(i-1,2) && x<=tar(i,2)
%         idx=i-1;
%         if x>=tar(i,1)
%             isInInt=1;
%         else
%             isInInt=0;
%         end
%         return
%     end
% end
% end
% 
% function [idx,isInInt]=getnextInt(tar,x)
% % idx = idx of next Int \in [1 size(tar,1)+1]
% % isInInt = 0 - no, 1 - yes with index idx-1
% if x>=tar(end,1)
%     idx=size(tar,1)+1;
%     if x<=tar(end,2)
%         isInInt=1;
%     else
%         isInInt=0;
%     end
%     return;
% end
% if x<tar(1,1)
%     idx=1;
%     isInInt=0;
%     return;
% end
% for i=1:size(tar,1)-1
%     if x>=tar(i,1) && x<tar(i+1,1)
%         idx=i+1;
%         if x<=tar(i,2)
%             isInInt=1;
%         else
%             isInInt=0;
%         end
%         return
%     end
% end
% end

