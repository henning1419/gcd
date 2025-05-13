function buttonTargetDown( hObject, eventdata,obj )
% set(obj.mainfig,'WindowButtonUpFcn',{@buttonTargetUp,obj},'WindowButtonMotionFcn',{@buttonTargetMove,obj},'KeyPressFcn','')
% 
% if ishandle(obj.activeplot)
%     delete(obj.activeplot)
% end
% 
% loc=get(obj.mainaxes, 'Currentpoint');
% obj.target.start=loc(1,1);
% obj.target.end=loc(1,1);
% 
% %initial plot
% yy=get(obj.mainaxes,'YLim');
% hold(obj.mainaxes,'on');
% obj.tarplot=fill(obj.mainaxes,obj.target.start*[1 1 1 1],...
%     [yy(1) yy(1) yy(2) yy(2)],...
%     obj.tarcol,...
%     'FaceAlpha',0.5);
% hold(obj.mainaxes,'off');

roi = drawrectangle(obj.mainaxes,'StripeColor','y')

end

