function [x, y, N] = buildIllum(patNum, N, dl, args)

% Effective "width" of radius, assume this is intrinsic beam sigma
radiusWidth = 0.025; 
totalNumPoints = 10000;


% estimated area fraction
estAreaFrac = 1;

% flip on raster
bFlipOnRaster = true;

switch patNum
    case 0 % annular/disk
        sigmaLow = args(1);
        sigmaHigh = args(2);
        
        mask = @(r, th) r > sigmaLow & r < sigmaHigh;
        bFlipOnRaster = false;
         
    case 1 % wedge pole 
        sigmaLow = args(1);
        sigmaHigh = args(2);
        wedgeWidth = args(3) * pi/180;
        mask = @(r, th) r > sigmaLow & r < sigmaHigh & th > -wedgeWidth/2 & th < wedgeWidth/2;
        
    case 2 % circular pole
        poleOffset = args(1);
        poleRadius = args(2);
        
        mask = @(r, th) (r.*cos(th) - poleOffset).^2 + (r.*sin(th)).^2 < poleRadius^2;
    
    case 3 % leaf dipole pole
        
        poleOffset = 2 - args(1);
        poleRadius = 1;
        mask = @(r, th) (r.*cos(th) - poleOffset).^2 + (r.*sin(th)).^2 < poleRadius^2;

 
    case 4 % leaf quasar pole
        deltaSigma =  args(1);
        poleRadius = 1;
        mask = @(r, th) (r.*cos(th) - (2 - deltaSigma)).^2 + (r.*sin(th)).^2 < poleRadius^2 & ...
           (r.*cos(th) - (1 - deltaSigma/2)).^2 + (r.*sin(th) - (1 - deltaSigma/2)).^2 < poleRadius^2 & ...
           (r.*cos(th) - (1 - deltaSigma/2)).^2 + (r.*sin(th) + (1 - deltaSigma/2)).^2 < poleRadius^2;
        
    case 5 % leaf hexapole pole
        deltaSigma =  args(1);
        poleRadius = 1;

        mask = @(r, th) ...
            (r.*cos(th + pi/6) - (2 - deltaSigma)).^2 + (r.*sin(th + pi/6)).^2 < poleRadius^2 & ...
            (r.*cos(th - pi/6) - (2 - deltaSigma)).^2 + (r.*sin(th - pi/6)).^2 < poleRadius^2;%& ...
%             (r.*cos(th - 0) - (2 - deltaSigma)).^2 + (r.*sin(th - 0) + (2 - deltaSigma)).^2 < poleRadius^2;
        
end

% Generate rdrdth grid:
radii = radiusWidth : radiusWidth : 1 - radiusWidth;
numRadii = length(radii);

% Estimate nominal dl, the point separation in theta.  This is constant for
% all rings to create equal density.  To see why this needs to be the case,
% consider that we want a constant dA = rdrdTh.  This means if dr is
% constant, dTh must be C/r, where C is a constant.  Alternatively, we
% require rdTh to be constant, which is dl here

if (dl == 0)
    dl = pi * (numRadii + 1)/totalNumPoints * 6 / estAreaFrac;
end



x = [];
y = [];
% concat all [X,Y]_r, flipping each row:

% loop through each radius, mask it, then concat
for  k = 1:length(radii)
    
    % Build a theta array for each radius
    dTh = dl/(2*pi*radii(k));
    thGrid = -pi : dTh : pi - dTh;
    
    % Mask this theta array to valid domain
    coordMsk = mask(radii(k) * ones(1, length(thGrid)) ,thGrid);
    [Xmsk, Ymsk] = pol2cart(thGrid(coordMsk(:)), radii(k));

    if ~isempty(Xmsk)
        
        if (mod(k,2) == 0 && bFlipOnRaster)
            x = [x; flipud(Xmsk(:))];
            y = [y; flipud(Ymsk(:))];
        else
            x = [x; Xmsk(:)];
            y = [y; Ymsk(:)];
        end
    end
    
end







