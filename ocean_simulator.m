function h = ocean_simulator(~)
% function h = ocean_simulator(~)
%
% Mini gui to experiment a few parameters used to create a wave simulation
% v1 - 2015/06/02 - Hoki 
% MIT licence
%
% This file is shared at the following link:
% https://www.dropbox.com/s/qow7cdf5z95t7hx/ocean_simulator.m?dl=0
%

%% // Default parameters
param.meshsize  = 128 ;     %// main grid size
param.patchsize = 200 ;     
param.windSpeed = 100  ;    %// what unit ? [m/s] ??
param.winddir   = 90   ;    %// Azimuth
param.rng = 13 ;            %// setting seed for random numbers
param.A         = 1e-7 ;    %// Scaling factor
param.g         = 9.81 ;    %// gravitational constant

param.xLim = [-10 10] ;     %// domain limits X
param.yLim = [-10 10] ;     %// domain limits Y
param.zLim = [-1e-4 1e-4]*2 ;

%% // Pre-set modes
presetModes.names     = {'Calm','Windy lake','Swell','T-Storm','Tsunami','Custom'} ;
presetModes.winddir   = [  90   ,    90   ,   220,     90     ,   90   ] ;
presetModes.patchsize = [  150  ,   500 ,    180   ,     128     ,   128   ] ;
presetModes.meshsize  = [  128  ,   480 ,    215   ,     128     ,   128   ] ;

%% // start first instance
if nargin == 1
    h  = init_gui(param,presetModes,true) ;	%// Create a new instance
else
    h  = init_gui(param,presetModes) ;	%// Delete all instances then redraw a new one
end
init_surf(h,param) ;                %// update the surface with default parameters

%% // save parameters into application data
setappdata(h.fig,'param',param) ; 
setappdata(h.fig,'presetModes',presetModes) ;


end

function [H0, W, Grid_Sign] =  initialize_wave( param )
% function [H0, W, Grid_Sign] =  initialize_wave( param )
%
% This function return the wave height coefficients H0 and W for the
% parameters given in input. These coefficients are constant for a given
% set of input parameters.
% Third output parameter is optional (easy to recalculate anyway)

    rng(param.rng);  %// setting seed for random numbers

    gridSize = param.meshsize * [1 1] ;

    meshLim = pi * param.meshsize / param.patchsize ;
    N = linspace(-meshLim , meshLim , param.meshsize ) ;
    M = linspace(-meshLim , meshLim , param.meshsize ) ;
    [Kx,Ky] = meshgrid(N,M) ;

    K = sqrt(Kx.^2 + Ky.^2);    %// ||K||
    W = sqrt(K .* param.g);     %// deep water frequencies (empirical parameter)

    [windx , windy] = pol2cart( deg2rad(param.winddir) , 1) ;

    P = phillips(Kx, Ky, [windx , windy], param.windSpeed, param.A, param.g) ;
    H0 = 1/sqrt(2) .* (randn(gridSize) + 1i .* randn(gridSize)) .* sqrt(P); % height field at time t = 0

    if nargout == 3
        Grid_Sign = signGrid( param.meshsize ) ;
    end

end

function Z = calc_wave( H0,W,time,Grid_Sign )
% Z = calc_wave( H0,W,time,Grid_Sign )
%
% This function calculate the wave height based on the wave coefficients H0
% and W, for a given "time". Default time=0 if not supplied.
% Fourth output parameter is optional (easy to recalculate anyway)

    % recalculate the grid sign if not supplied in input
    if nargin < 4
        Grid_Sign = signGrid( param.meshsize ) ;
    end
    % Assign time=0 if not specified in input
    if nargin < 3 ; time = 0 ; end
    
    wt = exp(1i .* W .* time ) ;
    Ht = H0 .* wt + conj(rot90(H0,2)) .* conj(wt) ;  
    Z = real( ifft2(Ht) .* Grid_Sign ) ;
end
    
function P = phillips(Kx, Ky, windDir, windSpeed, A, g)
%// The function now accept scalar, vector or full 2D grid matrix as input
    K_sq = Kx.^2 + Ky.^2;
    L = windSpeed.^2 ./ g;
    k_norm = sqrt(K_sq) ;
    WK = Kx./k_norm * windDir(1) + Ky./k_norm * windDir(2);
    P = A ./ K_sq.^2 .* exp(-1.0 ./ (K_sq * L^2)) .* WK.^2 ;
    P( K_sq==0 | WK<0 ) = 0 ;
end

function sgn = signGrid(n)
    [x,y] = meshgrid(1:n,1:n) ;
    sgn = ones( n ) ;
    sgn(mod(x+y,2)==0) = -1 ;
end

function init_surf(h,param)
%// this function update the surface graphic object with a new set of wave
%// parameter

%// Define the grid X-Y space
x = linspace( param.xLim(1) , param.xLim(2) , param.meshsize ) ;
y = linspace( param.yLim(1) , param.yLim(2) , param.meshsize ) ;
[X,Y] = meshgrid(x, y);

%// initialise wave coefficients
[H0, W, Grid_Sign] =  initialize_wave( param ) ;

%// calculate wave at t0
t0 = 0 ;
Z = calc_wave( H0 , W , t0 , Grid_Sign ) ;

%// Display the initial wave surface
set( h.surf , 'XData',X , 'YData',Y , 'ZData',Z )
set( h.ax   , 'XLim',param.xLim , 'YLim',param.yLim , 'ZLim',param.zLim )

%// Save coeffs for futur use
setappdata( h.fig , 'H0' , H0 )
setappdata( h.fig , 'W'  , W )
setappdata( h.fig , 'Grid_Sign'  , Grid_Sign )

setappdata( h.fig , 'param' , param) ;

%// Reset the wave animation parameters
animate_wave(h,true) ;

%// Display side patch if option checked
displaySidePatch = get(h.chkSidePatch,'Value') ; %// get checkbox state
init_SidePatch(h.fig,displaySidePatch) ; %// draw or erase side patches

end

function init_SidePatch(hobj,displaySidePatch)
    h = guidata( hobj ) ;
    param = getappdata( h.fig , 'param' ) ;

    if ~isfield( h , 'pt' ) ; h.pt = -1 ; end

    if displaySidePatch
        if any(~ishandle(h.pt))

            %// Delete the side patches in case some still exist
            delete(h.pt(ishandle(h.pt)))
%             if ishandle(h.sidewall) ; delete(h.sidewall) ; end
        
            %// create and display side patches
            x = linspace( param.xLim(1) , param.xLim(2) , param.meshsize ) ;
            y = linspace( param.yLim(1) , param.yLim(2) , param.meshsize ) ;
            pcol = [0.2857,1,0.7143] ;

            xface = [ x param.xLim(2) param.xLim(1) ] ;
            yface = [ y param.yLim(2) param.yLim(1) ] ;
            face0 = zeros(size(xface)) ;
            faceZ = [zeros(size(x)) param.zLim(1) param.zLim(1)] ;

            h.pt(1) = handle( patch( xface , face0+param.yLim(1) , faceZ ,pcol ) ) ;
            h.pt(3) = handle( patch( xface , face0+param.yLim(2) , faceZ ,pcol ) ) ;
            h.pt(2) = handle( patch( face0+param.xLim(1) , yface , faceZ ,pcol ) ) ;
            h.pt(4) = handle( patch( face0+param.xLim(2) , yface , faceZ ,pcol ) ) ;

            h.pt = handle(h.pt) ;
            set(h.pt, 'Facecolor',pcol , 'FaceAlpha',0.5 , 'EdgeColor','none')

%             xmin = param.xLim(1) ; xmax = param.xLim(2) ;
%             ymin = param.yLim(1) ; ymax = param.yLim(2) ;
%             zmin = param.zLim(1) ;
%             
%             lf = param.meshsize ; %facelength
%             lft = 4*lf ;
%             f0 = zeros(size(x)) ;
%             Xp = [x       f0+xmax fliplr(x) f0+xmin   xmin xmax xmax xmin ] ;
%             Yp = [f0+ymin y       f0+ymax   fliplr(y) ymin ymin ymax ymax ] ;
%             Zp = [f0      f0        f0        f0      zmin zmin zmin zmin ] ;
%             patchinfo.Vertices = [Xp.' , Yp.' , Zp.'] ;
%             patchinfo.Faces = [ ...
%                 1:lf        lft+2 lft+1 ; ...
%                 lf+1:2*lf   lft+3 lft+2 ; ...
%                 2*lf+1:3*lf lft+4 lft+3 ; ...
%                 3*lf+1:4*lf lft+1 lft+4] ;
%             patchinfo.FaceColor = pcol ;
% 
%             h.sidewall = handle( patch(patchinfo) ) ;   
%             set(h.sidewall, 'Facecolor',pcol , 'FaceAlpha',0.5 , 'EdgeColor','none')

            guidata( h.fig , h )
        end
        
        %// update top edge of patches
        Z = get( h.surf,'ZData') ;
        h.pt(1).ZData(1:end-2) = Z(1,:) ;
        h.pt(2).ZData(1:end-2) = Z(:,1) ;
        h.pt(3).ZData(1:end-2) = Z(end,:) ;
        h.pt(4).ZData(1:end-2) = Z(:,end) ;
%         ptop = [ Z(1,:).' ; Z(:,1) ; flipud(Z(end,:).') ; flipud(Z(:,end)) ; zmin ; zmin ; zmin ; zmin ] ;
%         h.sidewall.vertices(:,3) = ptop ;

    else %// = IF NOT displaySidePatch
        %// Delete the side patches
        delete(h.pt(ishandle(h.pt)))
%         if ishandle(h.sidewall) ; delete(h.sidewall) ; end
    end
    
end

function recalc_surf(hobj,~)
    h = guidata(hobj) ;
    param = getappdata(h.fig,'param') ;
    param = update_options( h , param ) ; %// Get settings parameters from sliders instead of hard coded
    init_surf( h , param)
end

function animate_wave(h,~)
    persistent time H0 W Grid_Sign hs
    
    if nargin == 2
    % Just re-set everything (has to be called once to initialise)
        time = 0 ;
        H0 = getappdata( h.fig , 'H0' ) ;
        W  = getappdata( h.fig , 'W'  ) ;
        Grid_Sign = getappdata( h.fig , 'Grid_Sign' ) ;      
        hs = handle( h.surf ) ;
        return
    end
    time = time + 0.1 ;
    
    %// update wave surface
    Z = calc_wave( H0,W,time,Grid_Sign ) ;
    hs.ZData = Z ;
    
    %// update side patches
    isPatchDisplayed = logical( get(h.chkSidePatch,'Value') ) ;
    if isPatchDisplayed
            if ishandle(h.pt(1)) ;  h.pt(1).ZData(1:end-2) = Z(1,:)   ;  end
            if ishandle(h.pt(2)) ;  h.pt(2).ZData(1:end-2) = Z(:,1)   ;  end
            if ishandle(h.pt(3)) ;  h.pt(3).ZData(1:end-2) = Z(end,:) ;  end
            if ishandle(h.pt(4)) ;  h.pt(4).ZData(1:end-2) = Z(:,end) ;  end
%     if ishandle(h.sidewall)
%         ptop = [ Z(1,:).' ; Z(:,1) ; flipud(Z(end,:).') ; flipud(Z(:,end)) ] ;
%         h.sidewall.vertices(1:end-4,3) = ptop ;
%     end
    end

end

function param = update_options( h , param )
%// read the values from the sliders
    param.meshsize  = round( get( h.sldMeshSize ,'Value' ) ) ;
    param.patchsize = round( get( h.sldPatchSize,'Value' ) ) ;
    param.winddir = get( h.sldWindDir,'Value' )  ;
    
    %// update the display
    set( h.lblMeshSize  , 'String', sprintf('Mesh size = %d'  ,param.meshsize )  ) ;
    set( h.lblPatchSize , 'String', sprintf('Patch size = %d' ,param.patchsize ) ) ;
    set( h.lblWindDir   , 'String', sprintf('Wind Dir = %g'   ,param.winddir )   ) ;
end

function h = init_gui(param,presetModes,~)

% Define a few dimensions reused accross the GUI (normalized units)
pnlH    = 0.2 ;  %// setting panel height
marginH = 0.1 ;  %// Horizontal margin
marginW = 0.05 ; %// Vertical margin
lineH   = 0.10 ; %// uicontrol line Height
lblW    = 0.2  ; %// label width
btnW    = 0.12 ; %// width for top line uicontrols

if nargin < 3
    hFigOld = findobj('Type','Figure','Tag','OceanSim') ;
    close(hFigOld)
end

h.fig   = figure('Name','Ocean simulator','Tag','OceanSim','CloseRequestFcn',@CloseGUI) ;
% set(h.fig,'WindowStyle','Docked')
% set(h.fig,'Position',[2074 248 1200 800])
movegui(h.fig,'center') 

h.timer = timer('Period',0.1,'TimerFcn',{@timerCallback,h.fig},'TasksToExecute',Inf ,'ExecutionMode','fixedRate') ;

h.pnlsettings = uipanel('Position',[0   0   1 pnlH   ]) ; %// lower panel for slider settings
h.pnldisplay  = uipanel('Position',[0 pnlH  1 1-pnlH ]) ; %// top panel for surface display

%% // populate the settings panel
%// common options for uicontrols
btnOptions = { 'Parent',h.pnlsettings , 'Style','pushbutton' , 'Unit','Normalized'  } ;
chkOptions = { 'Parent',h.pnlsettings , 'Style','checkbox'   , 'Unit','Normalized'  } ;
cmbOptions = { 'Parent',h.pnlsettings , 'Style','popupmenu'  , 'Unit','Normalized'  } ;
edtOptions = { 'Parent',h.pnlsettings , 'Style','edit'       , 'Unit','Normalized'  } ;
lblOptions = { 'Parent',h.pnlsettings , 'Style','text'       , 'Unit','Normalized' , 'String','' } ;
sldOptions = { 'Parent',h.pnlsettings , 'Style','Slider'     , 'Unit','Normalized' , 'Callback',@recalc_surf } ;

%// position calculations helper functions
getLblPos = @(ipos) [marginW        (ipos+1)*marginH+(ipos-1)*lineH lblW               lineH ] ;
getSldPos = @(ipos) [2*marginW+lblW (ipos+1)*marginH+(ipos-1)*lineH 1-(3*marginW+lblW) lineH ] ;
getBtnPos = @(ipos,jpos) [jpos*marginW+(jpos-1)*btnW (ipos+1)*marginH+(ipos-1)*lineH btnW lineH ] ;
ipos = 0 ; jpos = 0 ;

%// Mesh size settings
ipos = ipos+1 ;
h.lblMeshSize  = uicontrol( lblOptions{:} , 'Position',getLblPos(ipos) ) ;
h.sldMeshSize  = uicontrol( sldOptions{:} , 'Position',getSldPos(ipos) , 'Min',32 , 'Max',512 , 'Value',param.meshsize ) ;
%// Patch size settings
ipos = ipos+1 ;
h.lblPatchSize = uicontrol( lblOptions{:} , 'Position',getLblPos(ipos) ) ;
h.sldPatchSize = uicontrol( sldOptions{:} , 'Position',getSldPos(ipos) , 'Min',50 , 'Max',500 , 'Value',param.patchsize  ) ;
% %// Wind direction (Azimuth)
ipos = ipos+1 ;
h.lblWindDir  = uicontrol( lblOptions{:} , 'Position',getLblPos(ipos) ) ;
h.sldWindDir  = uicontrol( sldOptions{:} , 'Position',getSldPos(ipos) , 'Min',0 , 'Max',359 , 'Value',param.winddir , 'SliderStep',[1 10]/359 ) ;
%// Action buttons
ipos = ipos+1 ;
jpos = jpos+1 ; h.cmbPreset    = uicontrol( cmbOptions{:} , 'Position',getBtnPos(ipos,jpos) , 'String',presetModes.names , 'Value',1, 'Callback',@choose_preset ) ;
jpos = jpos+1 ; h.chkSidePatch = uicontrol( chkOptions{:} , 'Position',getBtnPos(ipos,jpos) , 'String','Side patch' , 'Value',0, 'Callback',@toggle_side_patch ) ;
jpos = jpos+1 ; h.btnAnimate   = uicontrol( btnOptions{:} , 'Position',getBtnPos(ipos,jpos) , 'String','Animate' , 'Callback',@toggle_animation ) ;
jpos = jpos+1 ; h.btnMakeGIF   = uicontrol( btnOptions{:} , 'Position',getBtnPos(ipos,jpos) , 'String','Make GIF' , 'Callback',@makeGIF ) ;
jpos = jpos+1 ; h.edtGIFname   = uicontrol( edtOptions{:} , 'Position',getBtnPos(ipos,jpos) , 'String','DancingWave.gif' ) ;


%// populate the display panel
h.ax   = axes('Parent',h.pnldisplay) ;  %// create the axes
h.surf = surf( NaN(2) ) ;               %// create an empty "surf" object
h.pt = -1 ;                             %// create an empty "patch" object placeholder
% h.sidewall = [] ;                     %// create an empty "patch" object placeholder (discontinued)

%// get actual object handles instead of only their pointer index
h.ax   = handle( h.ax ) ;
h.surf = handle( h.surf ) ; 

%// Define rendering options
h.ax.Position = get( h.pnldisplay,'Position') ;
axis off                                %// make the axis grid and border invisible
shading interp                          %// improve shading (remove "faceted" effect
%// create blue colormap
blue = linspace(0.4, 1.0, 25).' ; cmap = [blue*0, blue*0, blue];
colormap(cmap)
h.light_handle = lightangle(-45,30) ;   %// add a light source
%// configure lighting
set(h.surf,'FaceLighting','phong',...
    'AmbientStrength',.3,'DiffuseStrength',.8,...
    'SpecularStrength',.9,'SpecularExponent',25,...
    'BackFaceLighting','unlit')

guidata( h.fig , h ) % always useful to store your handle structure

end

function toggle_side_patch(hobj,~)
    displaySidePatch = get(hobj,'Value') ; %// get checkbox state
    init_SidePatch(hobj,displaySidePatch) ; %// draw or erase side patches
end

function toggle_animation(hobj,~)
    h = guidata( hobj ) ;
    btnStr = get(h.btnAnimate,'String') ; %// get button state
    if strcmp('Animate',btnStr)
        set( h.btnAnimate,'String','Stop' )
        start(h.timer) ;
    else
        set( h.btnAnimate,'String','Animate' )
        stop(h.timer) ;
    end
end

function choose_preset(hobj,~)
    h = guidata( hobj ) ;
    presetModes = getappdata(h.fig,'presetModes') ;

    psetIndex = round( get(h.cmbPreset,'Value') ) ;
    set(h.sldWindDir   ,'Value' , presetModes.winddir(psetIndex)) 
    set(h.sldPatchSize ,'Value' , presetModes.patchsize(psetIndex))
    set(h.sldMeshSize  ,'Value' , presetModes.meshsize(psetIndex))
    recalc_surf(h.fig) ;
    
end

function timerCallback(~,~,hfig)
    h = guidata( hfig ) ;
    animate_wave(h) ;
end

function CloseGUI(hfig,~)
    h = guidata(hfig) ;
    try delete(h.timer) ; catch ; disp('Timer could not be deleted') ; end
    delete(h.fig) % delete figure anyway
end

function makeGIF(hobj,~)
    h = guidata( hobj ) ;
    %// check for file extension
    gifname = get( h.edtGIFname,'String') ;
    if ~strcmpi( gifname(end-3:end) , '.gif' )
        gifname = [gifname '.gif'] ;
    end
    
    nFrame = 50 ;
    hframe = h.ax ;
    hsurf  = h.surf ;
    
%     time = 0 ;
    H0 = getappdata( h.fig , 'H0' ) ;
    W  = getappdata( h.fig , 'W'  ) ;
    Grid_Sign = getappdata( h.fig , 'Grid_Sign' ) ;      

    f = getframe(hframe) ;
    [im,map] = rgb2ind(f.cdata,256,'nodither');
    im(1,1,1,nFrame) = 0;
    iframe = 0 ;
    
    %// Alpha channel not supported in GIF files
    isPatchDisplayed = logical( get(h.chkSidePatch,'Value') ) ;
    if isPatchDisplayed
        alphapt = get( h.pt(ishandle(h.pt) ) , 'FaceAlpha') ;   %// save alpha value
        alphapt = alphapt{1} ;
        set( h.pt(ishandle(h.pt) ) , 'FaceAlpha',1 )            %// set to solid color
    end
    
    for time = (0:nFrame-1)*.5
        %// update wave surface
        Z = calc_wave( H0,W,time,Grid_Sign ) ;
        hsurf.ZData = Z ;

        %// update side patches
        if isfield( h , 'pt' )
            
            if ishandle(h.pt(1)) ;  h.pt(1).ZData(1:end-2) = Z(1,:)   ;  end
            if ishandle(h.pt(1)) ;  h.pt(2).ZData(1:end-2) = Z(:,1)   ;  end
            if ishandle(h.pt(1)) ;  h.pt(3).ZData(1:end-2) = Z(end,:) ;  end
            if ishandle(h.pt(1)) ;  h.pt(4).ZData(1:end-2) = Z(:,end) ;  end
        end
        pause(0.001);

        f = getframe(hframe) ;
        iframe= iframe+1 ;
        im(:,:,1,iframe) = rgb2ind(f.cdata,map,'nodither');
    end
    imwrite(im,map,gifname,'DelayTime',0,'LoopCount',inf) %g443800
    disp([num2str(nFrame) ' frames written in file: ' gifname])
    
    if isPatchDisplayed  %// restore alpha value for patch
        set( h.pt(ishandle(h.pt) ) , 'FaceAlpha',alphapt )            
    end

end