for i = 1:npoints
    xyvalues(i) = Settings.ImageNamesList(Settings.CalibrationPointIndecies(i));
end
xyvalues = xyvalues';
