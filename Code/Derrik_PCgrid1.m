% code to map out xstar, ystar, zstar for a load of points in the scan.
%keyboard PCMinSinglePattern and then run this (choose points near top left
%of IQ scan map to get close to 0,0 for calibration)
% maybe throw out, or weight less, points more than x-SDs from minimum point?
% Or based upon IQ of scan image? (I filter out values above a certain
% cutoff and then take median of the remaining points, and that works well
% on Area 1 steel)
tic
numpats=100;  % number of points to take from scan
numpc=40;    % number of PC points to take in each dimension (about 2 points per hour for a grid of 10x10x10 PCs)
deltapc=0.06/numpc;   %step size for PC grid search
ScanParams.xstar=0.5108; %local minima for BYU Martensite steel
ScanParams.ystar=0.5963;
ScanParams.zstar=0.6702;
q=input('x=1,y=2 or z=3 variable for PC ');
 %counter = 0;
    
for qq=1:numpats
    qq
    Ind=round(Settings.ScanLength/numpats*qq);
    Indsave(qq)=Ind;
    %Apply Plane Fit
    ConfIndex(qq) = Settings.CI(Ind);
    Fitplot(qq) = Settings.Fit(Ind);
    
    %Try to filter out useless points based on CI and Fit
    if ConfIndex(qq) > .6 && Fitplot(qq) < .8
         %counter = counter +1;
    
    xstar = ScanParams.xstar-Settings.XData(Ind)/Settings.PhosphorSize;
    ystar = ScanParams.ystar+Settings.YData(Ind)/Settings.PhosphorSize*sin(Settings.SampleTilt-Settings.CameraElevation);
    zstar = ScanParams.zstar+Settings.YData(Ind)/Settings.PhosphorSize*cos(Settings.SampleTilt-Settings.CameraElevation);
    PC0(1) = xstar;
PC0(2) = ystar;
PC0(3) = zstar;
%Set up a CI vector like ci = settings.CI(Ind)
if q==1
    star=xstar;
elseif q==2
    star=ystar;
else
    star=zstar;
end
    Av = Settings.AccelVoltage*1000; %put it in eV from KeV
    sampletilt = Settings.SampleTilt;
    elevang = Settings.CameraElevation;
    pixsize = Settings.PixelSize;
    Material = ReadMaterial(Settings.Phase{Ind});
    
    % keyboard
    ImagePath = Settings.ImageNamesList{Ind};
    ScanImage = ReadEBSDImage(ImagePath,Settings.ImageFilter);
    
    [roixc,roiyc]= GetROIs(ScanImage,Settings.NumROIs,pixsize,Settings.ROISize,...
        Settings.ROIStyle);
    Settings.roixc = roixc;
    Settings.roiyc = roiyc;
    g = euler2gmat(Settings.Angles(Ind,1),Settings.Angles(Ind,2),Settings.Angles(Ind,3)); % DTF - don't use ref angles for grain as is done on previous line!!
    
    for xx=1:numpc
                PC0(q) = star+(xx-1-(numpc-1)/2)*deltapc;
               paramspat={PC0(1);PC0(2);PC0(3);pixsize;Av;sampletilt;elevang;Material.Fhkl;Material.dhkl;Material.hkl};
                PCvals(xx)=PC0(q);
                pctest(qq,xx)=CalcNormFMod(PC0,ScanImage,paramspat,Material.lattice,Material.a1,Material.b1,Material.c1,Material.axs,g,Settings.ImageFilter,Ind,Settings);
    end
    end
end

toc


figure;hold on ;
for i=1:numpats; plot(PCvals,squeeze(pctest(i,:)),'*');end

%find out how many values are 'reasonable'
clear lowpctest
for i=1:numpc;numlow(i)=sum(pctest(:,i)<5e-3);end
[nlow indlow]=max(numlow);
% pick out the lowest values for each pc point
for i=1:numpc
    temp=sort(pctest(:,i));
    lowpctest(:,i)=temp(1:nlow);
end
figure;hold on ;
for i=1:nlow; plot(PCvals,squeeze(lowpctest(i,:)),'*');end   
for i=1:numpc; thismean(i)=median(lowpctest(:,i))    ;end
plot(PCvals,thismean);ylim([0 .001]);

% alternatively - pick out the patterns corresponding to those with low
% strain values at on pc
qq=(pctest(:,indlow)<5e-3);%Could change this value as well if needed
nn=[1:numpats];
nn=nn(qq);
figure;hold on ;
for i=1:nlow; plot(PCvals,(pctest(nn(i),:)),'*'); end   
for i=1:numpc; thismean(i)=median(pctest(qq,i))    ;end
plot(PCvals,thismean); ylim([0 .01]);

% fit a quadratic to the data to find minimum
pp=polyfit(PCvals,thismean,2);
PCopt=-pp(2)/2/pp(1) % find minimum for optimal PCx (note that this is for Given PCy and PCz, so may not be the optimal for all 3)


return
figure;hold on ;
for i=1:numpats; plot(PCyvals,squeeze(pctest(i,7,:,3)),'*');end
for i=1:numpc; pctemp=squeeze(pctest(:,7,i,3)); thismean(i)=mean(pctemp(pctemp<2e-3))    ;end
for i=1:numpc; pctemp=squeeze(pctest(:,7,i,3)); thismean(i)=median(pctemp(pctemp<3e-3))    ;end
plot(PCyvals,thismean);
pp=polyfit(PCyvals,thismean,2);
PCyopt=-pp(2)/2/pp(1)

figure;hold on ;
for i=1:numpats; plot(PCzvals,squeeze(pctest(i,7,10,:)),'*');end
for i=1:numpc; pctemp=squeeze(pctest(:,7,10,i)); thismean(i)=mean(pctemp(pctemp<2e-3))    ;end
for i=1:numpc; pctemp=squeeze(pctest(:,7,10,i)); thismean(i)=median(pctemp(pctemp<3e-3))    ;end
plot(PCzvals,thismean);
pp=polyfit(PCzvals,thismean,2);
PCzopt=-pp(2)/2/pp(1)

% For area 1 steel:
% PCxopt =   0.521642746394567
% PCyopt =   0.729096609966661
% PCzopt =   0.528060416859464

% to plot surface of median values
for i=1:numpc 
    for j=1:numpc
    pctemp=squeeze(pctest(:,i,j,3)); 
    surfmed(i,j)=median(pctemp(pctemp<2e-3))    ;
    end
end
figure
surf(surfmed)


% now run Fminsearch at all the points and see if the average results
% agrees with the minimum of the CalcNormFMod surface:

PCarray=zeros(numpats,3);
for qq=1:numpats
    qq
    Ind=round(Settings.ScanLength/numpats*qq);
    Indsave(qq)=Ind;
    %Apply Plane Fit
    xstar = ScanParams.xstar-Settings.XData(Ind)/Settings.PhosphorSize;
    ystar = ScanParams.ystar+Settings.YData(Ind)/Settings.PhosphorSize*sin(Settings.SampleTilt-Settings.CameraElevation);
    zstar = ScanParams.zstar+Settings.YData(Ind)/Settings.PhosphorSize*cos(Settings.SampleTilt-Settings.CameraElevation);
    
    PC0(1) = xstar;
    PC0(2) = ystar;
    PC0(3) = zstar;
    
    Av = Settings.AccelVoltage*1000; %put it in eV from KeV
    sampletilt = Settings.SampleTilt;
    elevang = Settings.CameraElevation;
    pixsize = Settings.PixelSize;
    Material = ReadMaterial(Settings.Phase{Ind});
    
    % keyboard
    ImagePath = Settings.ImageNamesList{Ind};
    ScanImage = ReadEBSDImage(ImagePath,Settings.ImageFilter);
    
    [roixc,roiyc]= GetROIs(ScanImage,Settings.NumROIs,pixsize,Settings.ROISize,...
        Settings.ROIStyle);
    Settings.roixc = roixc;
    Settings.roiyc = roiyc;
    g = euler2gmat(Settings.Angles(Ind,1),Settings.Angles(Ind,2),Settings.Angles(Ind,3)); % DTF - don't use ref angles for grain as is done on previous line!!
    options = optimset('TolX',1e-6,'TolFun',1e-6);
    [PCprime,value,flag,iter] = fminsearch(@(PC)CalcNormFMod(PC,ScanImage,paramspat,Material.lattice,Material.a1,Material.b1,Material.c1,Material.axs,g,Settings.ImageFilter,Ind,Settings),PC0,options);
    PCarray(qq,:)=PCprime(:);
end

figure
plot(PCarray(:,1),PCarray(:,2),'*')
hold on
plot(mean(PCarray(:,1)),mean(PCarray(:,2)),'r*')
% area 1 PCx and y values: 0.523429726772288, 0.726952852711562


% for jj = 1:100
%     ConfIndex(jj) = Settings.CI(Ind);
% end
% %ConfIndex(qq) = Settings.CI(Ind);
% figure(1)
% plot(ConfIndex, pctest(:,1),'*')