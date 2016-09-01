xdata = Settings.CalibrationPointsPC(:,1)+(Settings.XData(Settings.CalibrationPointIndecies))/psize;
ydata = Settings.CalibrationPointsPC(:,2)-(Settings.YData(Settings.CalibrationPointIndecies))/psize*sin(Settings.SampleTilt);
zdata = Settings.CalibrationPointsPC(:,3)-(Settings.YData(Settings.CalibrationPointIndecies))/psize*cos(Settings.SampleTilt);

meanx = mean(xdata);
meany = mean(ydata);
meanz = mean(zdata);

figure(1)
plot(xdata,ydata,'bo')
hold on
plot(meanx,meany,'rx')