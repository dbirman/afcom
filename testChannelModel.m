%% Test channel model

parfor i = 1:100
    %% Build 360 units, 180 left tuned, 180 right tuned

    n = 1000;
    n2 = n/2;
    % each unit is represented by three columns, the first column is the
    % spatial tuning [-1,1], the second and third column represented the center
    % of a von mises distribution corresponding to color/motion sensitivity:
    units = zeros(n,3);
    % set spatial
    units(1:n2,1) = -1;
    units((n2+1):end,1) = 1;
    % set color sensitivity
    units(1:n2,2) = rand(n2,1)*2*pi;
    units((n2+1):end,2) = rand(n2,1)*2*pi;
    % set motion sensitivity
    units(1:n2,3) = rand(n2,1)*2*pi;
    units((n2+1):end,3) = rand(n2,1)*2*pi;
    % for plotting, compute the x/y of color and motion sensitivity
    xc = cos(units(:,2));
    yc = sin(units(:,2));
    xm = cos(units(:,3));
    ym = sin(units(:,3));

    %% Show design matrix
    figure;
    imagesc(units);

    %% Parameters
    cK = 1; % kappa parameter for color
    mK = 1; % kappa parameter for motion

    %% Test a stimulus

    stim = [-1 pi pi]; 

    out = compTestChannelOut(units,stim,cK,mK);

    % plot two figures in color and motion space
    figure
    subplot(211)
    plot(out.*xc,out.*yc,'o','MarkerFaceColor','k','MarkerEdgeColor','w');
    hline(0,'--k');
    vline(0,'--k');
    axis equal
    title('Color');
    subplot(212)
    plot(out.*xm,out.*ym,'o','MarkerFaceColor','k','MarkerEdgeColor','w');
    hline(0,'--k');
    vline(0,'--k');
    axis equal
    title('Motion');


    %% Create a four channel stimulus

    % up/down motion
    stims = [-1 rand*2*pi 0
             0 rand*2*pi pi
             0  rand*2*pi 0
             0  rand*2*pi pi];

    out = zeros(size(stims,1),size(units,1));
    for si = 1:size(stims,1)
        out(si,:) = compTestChannelOut(units,stims(si,:),cK,mK);
    end

    aout = sum(out,1)';

    figure
    subplot(211)
    plot(aout.*xc,aout.*yc,'o','MarkerFaceColor','k','MarkerEdgeColor','w');
    hline(0,'--k');
    vline(0,'--k');
    axis equal
    subplot(212)
    plot(aout.*xm,aout.*ym,'o','MarkerFaceColor','k','MarkerEdgeColor','w');
    hline(0,'--k');
    vline(0,'--k');
    axis equal

    %% Create a readout filter

    filt = [-1 nan 0]; % for example: left side and upward motion, what was the color
    fcK = 15;
    fmK = 15;

    %% Use the filter to estimate your probability of making different responses

    % multiply by the PDF of this filter 
    fout = compTestChannelOut(units,filt,fcK,fmK);

    figure;
    plot(fout);
    % fout = fout > 0.5;

    readout = aout .* fout;

    figure
    subplot(211)
    plot(readout.*xc,readout.*yc,'o','MarkerFaceColor','k','MarkerEdgeColor','k');
    hline(0,'--k');
    vline(0,'--k');
    axis equal
    subplot(212)
    plot(readout.*xm,readout.*ym,'o','MarkerFaceColor','k','MarkerEdgeColor','k');
    hline(0,'--k');
    vline(0,'--k');
    axis equal


    %% Compute estimate of location
    % these are actually angles originally, so go to X/Y space and average
    xread = readout .* xc;
    yread = readout .* yc;
    readang = atan2(mean(yread),mean(xread));
    disp(sprintf('True angle %1.2f rad',stims(1,2)));
    disp(sprintf('Readout angle %1.2f rad',readang));
    d(i) = angdist(stims(1,2),readang);
end