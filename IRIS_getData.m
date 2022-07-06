
function [Data]=IRIS_getData()

    FileName = uigetfile('*.xls','Selezionare file');
    fid = fopen(FileName);    
    n_line = 1;
    while 1==1                
    tline=fgetl(fid);            
        if tline==-1
           break
        end
        test=isempty(tline);
        if test==1
        elseif tline(1)=='A'
        
        else
            Data(:,n_line)=sscanf(tline,'%g',inf);
            n_line = n_line+1;
        end
    end


    fclose('all');
end

