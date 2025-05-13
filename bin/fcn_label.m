function fcn_label(hobj,event,obj,id)

disp('Result data for the clicked point:')
disp(['=> id:             ' num2str(id)])
disp(['=> profile length: ' num2str(obj.ds_length(id))])
disp(['=> peak count:     ' num2str(obj.ds_nr_peaks(id))])
disp(' ')

end