function usObj = readPXI5752Xml(usObj,file)

    [filepath,name,ext] = fileparts(file);

    fid = fopen(file);
    f = fread(fid,'*char')';
    fclose(fid);

    f = regexprep(f,'á','a');
    f = regexprep(f,'ê','e');
    fid = fopen(file,'w');
    fprintf(fid,'%s',f);
    fclose(fid);

    mystruct = util.xml2struct2(file);

    % initiate a object
    usObj = sensor.Ultrasonic;
    usObj.filename = name;
    usObj.filepath = filepath;
    usObj.fs = 0;
    usObj.fprf = 0;
    usObj.nwaves = 0;

    for i = 1 : size(mystruct.LVData.I32, 2)
        a = cell2mat(mystruct.LVData.I32(i));

        if strfind(a.Name.Text, 'Number')
            usObj.nwaves = str2num(a.Val.Text);
        elseif strfind(a.Name.Text, 'Pulse')
            usObj.fprf = str2num(a.Val.Text);  
        elseif strfind(a.Name.Text, 'Frequencia')
            usObj.fs = str2num(a.Val.Text);  
        elseif strfind(a.Name.Text, 'sampling')
            usObj.fs = str2num(a.Val.Text);  
        end  
    end
end

