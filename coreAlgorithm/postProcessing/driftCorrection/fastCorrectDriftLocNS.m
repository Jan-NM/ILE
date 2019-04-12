function OrteCorr = fastCorrectDriftLocNS(Orte,steps,volume,dbg)
% Corrects the drift of localization data where no permanent structure in
% the raw image is apparent and doesn't reveal permanent information about
% the sample's location.
%
% Synopsis: OrteCorr=fastCorrectDriftLocNS(Orte,steps,volume,dbg)
%
% Input:
%  Orte:        The 10xN Orte Matrix where N denotes the number of found
%               fluorophores
%  steps:       Number of steps to calculate the drift
%               default: 60
%  volume:      Number of slices considered to generate the images that are
%               correlated
%               Default: Value to cover all present fluorophes. Set to -1
%               for default.
%  dbg:         Debug. If true the images that are to be correlated and the
%               cross-correlation are displayed.
%               Default: false
% Output:
%  OrteCorr:    Drift corrected Orte Matrix.
%
% Drift correction for localization data only. Images are generated based
% on the Orte Matrix and the shift between these and a mean image taken from
% a subset located around the center of time-stack is determined using
% driftfindshift. 



%A 10th order spline is then fitted through the
% acquired shift estimations. The detected fluorphores provided by the Orte
% Matrix are then shifted based on this polynomial by simple vector-transformations.
% Iteration may lead to better results:
% OrteCorr = fastCorrectDriftLocNS(fastCorrectDriftLocNS(Orte,100),100)
%

%sampPoints and shiftvec - output


% Ideas based on findDrift by Rainer Kaufmann.
%
% Kirti Prakash
% k.prakash@imb-mainz.de
% Martin Hagmann
% martin.hagmann@kip.uni-heidelberg.de
% April 2014, Cremer Group
%
% modified by Jan Neumann - 26.09.2017
% added spline interpolation
% added RCC method
% changed initialization routine (e.g. calculation of stacks)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global SPDMparamstruct;

if nargin < 1
    error('Where do you want to perform the correction on? Provide at least an Orte Matrix.')
end

if nargin < 2
    steps = 30;
    warning('No number of samples provided. Defaulting to 60.')
end

if nargin < 4
    dbg = false;
end

if ~issorted(Orte(:,9))
    Orte = sortrows(Orte,9);
    warning('Detected that the Orte you provided is not sorted along the frame number. Did that for you.')
end

if nargin < 3 || volume < 0
    % Compute the sub stack size (volume) of one sample point so that all
    % volumes comprise the entire raw stack
    volume = floor((Orte(end,9)-Orte(1,9))/steps);
    %fprintf('Defaulted volume to %g frames to cover the entire orte matrix\n',volume)
elseif mod(volume,1)
    error('You need to specify the volume as an positive integer')
end

if ~isa(dbg,'logical')
    error('The debug value needs to be boolean.')
end

if steps < 2
    error('Need at least 2 steps to compute the drift.')
end

if steps >= length(unique(Orte(:,9)))
    steps = length(unique(Orte(:,9)))-1;
    warning(['There were more steps specified than frames in the Orte Matrix are present. Setting to ' num2str(steps) ', which is the amount of frames in your Orte Matrix.'])
end

fmax = max(Orte(:,9));
fmin = min(Orte(:,9));

% define parameters
pxlsz = 10;
th = 100; %maximum acceptable relative shift(x)+shift(y) in nm

%generate images at locations stored in sampPoints
sampPoints = round((Orte(1, 9):volume:Orte(end, 9)));
sampPoints = [fmin sampPoints(1, 1:end-1)+ceil(volume/2)  fmax];
%% find the index vectors (find corresponding indices in the Orte Matrix)
%allocation of the index vectors
steps = size(sampPoints, 2); % adapt steps to number of samp Points
indxFirst = zeros([steps 1]);
indxLast = indxFirst;

%if no entry is found at sampPoins +/- volume*.5, search for the
%next/previous frame, bring error if runs beyond the orte file
for i=2:steps-1 %first and last stacks are a bit different
    ispresent = false;
    counter = 0;
    while ~ispresent
        try %if a fluorophore is inside
            indxFirst(i)  = find(Orte(:,9)==sampPoints(i)-ceil(volume/2)-counter,1,'first');
            ispresent = true;
        catch %increment counter if not
            counter = counter+1;
            if sampPoints(i)-ceil(volume/2)-counter<Orte(1,9)
                error('Did not find sufficient entries in the Orte file. Try reducing the volume.')
            end
        end
    end
    ispresent = false;
    counter = 0;
    while ~ispresent
        try %if a fluorophore is inside
            indxLast(i) = find(Orte(:,9)==sampPoints(i)+ceil(volume/2)+counter,1,'last');
            ispresent = true;
        catch
            counter = counter+1;
            if sampPoints(i)+ceil(volume/2)+counter>Orte(end,9)
                error('Did not find sufficient entries in the Orte file. Try reducing the volume.')
            end
        end
    end
end

% determine the stacks for the first and last volume
% first stack is at first frame - drift is zero at this point
indxFirst(1) = find(Orte(:,9)==sampPoints(1),1,'first');
indxLast(1) = indxFirst(1);
% last stack is at last frame - drift is equal to the value of previous frame
indxLast(end) =  find(Orte(:,9)==sampPoints(end),1,'last');
indxFirst(end) = indxLast(end);

%% allocate the shiftvectors
shiftvec = zeros(steps, 2);
shiftvec(1,:) = [0 0];%MATLAB sucks. The first one means the zeroth which of course is zero

% waitbar
wait = waitbar(0,'determining drift...');
ph = findobj(wait, 'type', 'patch');
set(ph, 'FaceColor', [0.3 0.5 1], 'EdgeColor', [0 0 0]);
s=clock;
esttime = (pi+pi+pi+pi)^pi/(pi+pi);%we've been eating spacecake

%% do drift correction
% reference size image to make correlation easier
x1 = round(Orte(:,3)./pxlsz);
y1 = round(Orte(:,2)./pxlsz);

minX = min(x1);
minY = min(y1);
if max(y1-min(y1)) > max(x1-min(x1))
    [xSize, ySize] = size(full(sparse(y1-min(y1)+1,y1-min(y1)+1,1)));
else
   [xSize, ySize] = size(full(sparse(x1-min(x1)+1,x1-min(x1)+1,1)));
end
% direct cross-correlation
OrteRef = Orte(indxFirst(2):indxLast(2), :); % works only if the orte is sorted along the nineth column
imorigin = histogramBinnedDrift(OrteRef, pxlsz, xSize, ySize, minX, minY);
autocorr = driftfindshift(imorigin, imorigin);
[xCenter, yCenter] = findCorrelationPeak(autocorr);

for step = 3:steps-1
    subOrte = Orte(indxFirst(step):indxLast(step),:);
    imbim = histogramBinnedDrift(subOrte, pxlsz,  xSize, ySize, minX, minY);
    corr = driftfindshift(imorigin, imbim);
    [xShift, yShift] = findCorrelationPeak(corr);
    shiftvec(step, 1) = (xShift-xCenter)*pxlsz;
    shiftvec(step, 2) = (yShift-yCenter)*pxlsz;
    
    % threshold warning
    if sum(abs(shiftvec(step,:)-shiftvec(step-1,:)) > th)
        warning(['Found relative shift of [',num2str(shiftvec(step,:)-shiftvec(step-1,:)),']. This exceeds threshold of ', num2str(th) , ' nm in iteration step ', num2str(step), '.'])
    end
    
    if step==2
        s = clock;
    elseif step==3
        is = etime(clock,s);
        esttime = is*(steps-4);
    end
    wait = waitbar(step/(steps-1),wait,['remaining time: ',num2str(floor((esttime-etime(clock,s))/60), '%i'), ' minutes and ',num2str(floor(mod(esttime-etime(clock,s),60)), '%i'),' sec ...' ]);
end
shiftvec(2, :) = [0 0]; % first subset has zero shift
shiftvec(1, :) = shiftvec(2, :); % first frame has shiftvector equal to first subset
shiftvec(end, :) = shiftvec(end-1, :); % last frame has shiftvector equal to last subset
% mean shift
meanShift = sqrt((shiftvec(1, 1)-shiftvec(:, 1)).^2+(shiftvec(1, 2)-shiftvec(:, 2)).^2);


waitbar(100 ,wait, 'correcting shift in Orte Matrix. This may take a few minutes...');
%% fitting/interpolation and plotting of drift
% interpolate drift vector
xdrift = interp1(sampPoints, transpose(shiftvec(:,1)), fmin:fmax, 'spline');
ydrift = interp1(sampPoints, transpose(shiftvec(:,2)), fmin:fmax, 'spline');
meanDrift = interp1(sampPoints, transpose(meanShift(:, 1)), fmin:fmax, 'spline');
if SPDMparamstruct.driftCorrPlot == 1
    figure;
    hold on
    xlim([sampPoints(1) sampPoints(end)+10]);
    plot((sampPoints),shiftvec,'x');
    plot((fmin:fmax), transpose(xdrift));
    plot((fmin:fmax), transpose(ydrift));
    legend('x','y','Location','NorthEast');
    plot((fmin:fmax), transpose(meanDrift),'DisplayName','mean drift');
    xlabel('Frame number');
    ylabel('Drift in nm');
    % line([fmin, fmax], [95, 95]);
    title('Drift Correction Plot');
    hold off
end
%% save shift parameter
%sigx=sqrt(var(shiftvec(:,1)-transpose(polyval(px,sampPoints))));
%sigy=sqrt(var(shiftvec(:,2)-transpose(polyval(py,sampPoints))));
%save('shiftcorrparams.mat','px','py','sampPoints','shiftvec','sigx','sigy');
%% Correct the drift based on the obtained drift
%disp('correcting shift in the Orte Matrix...')
OrteCorr = Orte;%allocate

for ii = fmin:fmax
    idx = find(Orte(:, 9) == ii);
    OrteCorr(idx, 2) = Orte(idx, 2) + ydrift(1, (ii-fmin+1) ); % -fmin+1 to correct offset amongst ydrift and Orte, if Orte starts from frame number > 1
    OrteCorr(idx, 3) = Orte(idx, 3) + xdrift(1, (ii-fmin+1 ));
end
%% Repair the orte matrix because we might have got some negative locations
minY=min(OrteCorr(:,2));%y
minX=min(OrteCorr(:,3));%x
if minY<0
    OrteCorr(:,2)=OrteCorr(:,2)-minY;
end
if minX<0
    OrteCorr(:,3)=OrteCorr(:,3)-minX;
end
close(wait);
end