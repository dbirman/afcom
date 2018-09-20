%% Plot the L*a*b color space at certain luminance values by converting into RGB

msc = initScreen;

%% Plot the entire L*a*b color space at L==0
mglClearScreen(0);

as = -128:10:127;
bs = -128:10:127;

L = 53.3889647411143; % luminance value for gray screen

usewidth = msc.imageWidth/2;
useheight = msc.imageHeight/2;

% compute box widths
xwid = usewidth/length(as);
ywid = useheight/length(bs);

disppercent(-1/length(as));
for ai = 1:length(as)
    a = as(ai);
    for bi = 1:length(bs)
        b = bs(bi);
        
        rgb = lab2rgb([L a b],'ColorSpace','adobe-rgb-1998');
        
        degx = ai/length(as)*usewidth-usewidth/2;
        degy = bi/length(bs)*useheight-useheight/2;
        
%         disp(sprintf('Plotting L*a*b %1.1f %1.1f %1.1f: %1.0f %1.0f',L,a,b,degx,degy));
        
%         if
        mglFillRect(degx,degy,[xwid ywid],rgb);
%         end
    end
    disppercent(ai/length(as));
end
disppercent(inf);

mglFlush

%% Now plot the circle of values that are D distance away from the grey point
mglClearScreen([0.5 0.5 0.5]);

Lab = rgb2lab([0.5 0.5 0.5]);

L = Lab(1);
acenter = Lab(2);
bcenter = Lab(3);

D = 70;

usewidth = msc.imageWidth/2;
useheight = msc.imageHeight/2;

% compute box widths
xwid = usewidth/length(as);
ywid = useheight/length(bs);

% actual range of values used in experiment (note that the luminance seems
% off in sRGB mode... use adobe-rgb? 
theta_ = pi/64;
thetas = 0:theta_:2*pi;
dispthetas=thetas;
% thetas = [0 + (-pi/4:theta_:pi/4) pi + (-pi/4:theta_:pi/4)]; %0:theta_:2*pi;
% dispthetas = [0 + (-pi/2:theta_*2:pi/2) pi + (-pi/2:theta_*2:pi/2)];

offset = 0;%rand*2*pi;

for ti = 1:length(thetas)
    theta = thetas(ti);
    dtheta = offset+dispthetas(ti);
    % rotate through the Lab space plotting the colors at this Luminance
    % value
    
    a = D*cos(theta)+acenter;
    b = D*sin(theta)+bcenter;
    
    rgb = lab2rgb([L a b],'ColorSpace','adobe-rgb-1998');
    
    degx = ti/length(thetas)*usewidth-usewidth/2;
    
    mglGluPartialDisk(0,0,7.5,12.5,180/pi*(dtheta-theta_),180/pi*theta_*2,rgb);
%     mglFillRect(degx,0,[xwid xwid],rgb);
end

% mglGluPartialDisk(0,0,7.5,12.5,180/pi*(offset+pi/2-pi/64),180/pi*pi/32,[0.5 0.5 0.5]);
% mglGluPartialDisk(0,0,7.5,12.5,180/pi*(offset-pi/2-pi/64),180/pi*pi/32,[0.5 0.5 0.5]);


mglFlush