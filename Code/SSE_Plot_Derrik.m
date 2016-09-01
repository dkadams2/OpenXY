%Set up variables for plotting
CI = Settings.CI;
%  SSE1 = Settings.SSE;
%   SSE = cell2mat(SSE1);
Nx = Settings.Nx;
Ny = Settings.Ny;
ScanType = Settings.ScanType;

%Create Plots
        CIPlot = reshape(CI, Nx,Ny)';
        if Ny == 1 %Lines Scans
            CIPlot = repmat(CIPlot,floor(Settings.ScanLength/4),1);
        end

   
        imagesc(CIPlot)
        axis image %scales to natural width and height
%caxis([0 1]);
%colormap jet