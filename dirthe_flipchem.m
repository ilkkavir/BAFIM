function comp = dirthe_flipchem(Ne,Te,Ti,glat,glon,alt,fc)
%
% comp = dirthe_flipchem(Ne,Te,Ti,modErr,date,glat,glon,alt)
%
% ion composition from flipchem. modErr is error term in flipchem composition
% (0 to get the true model values)
%
%
% IV 2021

    % date in python format
%    pydate = py.datetime.datetime(int32(date.Year),int32(date.Month),int32(date.Day),int32(date.Hour),int32(date.Minute));

    % initialize flipchem
%    fc = py.flipchem.Flipchem(pydate);
    
    % call the model
    outputs = fc.get_point(glat,glon,alt,Ne,Te,Ti);

    % % convert the output list into a vector
    % outputsm = cell(outputs);

    % the final composition
    comp = outputs{4}/Ne;% + modErr;
    
end

