function buttonTargetMove( hObject, eventdata ,obj)

loc=get(obj.mainaxes, 'Currentpoint');
cloc=loc(1,1);
xx=get(obj.mainaxes,'XLim');
yy=get(obj.mainaxes,'YLim');

if isfield(obj.target,'start')
    sloc=obj.target.start;
    if cloc>xx(2)
        cloc=xx(2);
    end
    
    if cloc<xx(1)
        cloc=xx(1);
    end
    obj.target.end=cloc;

    set(obj.tarplot,'Xdata',[sloc cloc cloc sloc],'Ydata',[yy(1) yy(1) yy(2) yy(2)])

end
end

