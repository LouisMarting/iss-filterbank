function [h] = niceplot(varargin)
    
    
    % Parameters
    fontSize       = 18; %25
    lineWidthFrame = 1; 2;
    lineWidth      = 2; %3
    markerSize     = 3;
%     markerEdgeColor= 'none';
    fontName       = 'Times';
    
    extraArg=0;
    
    %% If there is no specification of lineWidth, use default of lineWidth
    lineWidth_flag = 0;
    for i=1:nargin
        if strcmpi('lineWidth',varargin{i})
            lineWidth_flag = 1;
        end
    end
    if ~lineWidth_flag
        varargin{nargin+extraArg+1}='lineWidth';
        varargin{nargin+extraArg+2}=lineWidth;
        extraArg=extraArg+2;
    end
    
    %% If there is no specification of markerSize, use default of markerSize
    markerSize_flag = 0;
    for i=1:nargin
        if strcmpi('markerSize',varargin{i})
            markerSize_flag = 1;
        end
    end
    if ~markerSize_flag
        varargin{nargin+extraArg+1}='markerSize';
        varargin{nargin+extraArg+2}=markerSize;
        extraArg=extraArg+2;
    end
    
    
     %% If there is no specification of markerEdgeColor, use default of markerEdgeColor
%     markerEdgeColor_flag = 0;
%     for i=1:nargin
%         if strcmpi('markerEdgeColor',varargin{i})
%             markerEdgeColor_flag = 1;
%         end
%     end
%     if ~markerEdgeColor_flag
%         varargin{nargin+extraArg+1}='markerEdgeColor';
%         varargin{nargin+extraArg+2}=markerEdgeColor;
%         extraArg=extraArg+2;
%     end
%     
    
    %% Plot
    h=plot(varargin{:});
    drawnow;
    
    %% Make nicer
    box on;
    set(gca,'XGrid','on');
    set(gca,'YGrid','on');
    set(gca,'ZGrid','on');
%     set(gca,'XMinorGrid','on');
%     set(gca,'YMinorGrid','on');
%     set(gca,'ZMinorGrid','on');
    set(gca,'GridLineStyle','-'); %Default '-'
    set(gca,'MinorGridLineStyle',':'); %Default ':'
    set(gca,'MinorGridColor',0.1*[1,1,1]); %Default 0.1
    set(gca,'GridColor',0.15*[1,1,1]);  %Default 0.15
    set(gca,'fontSize',fontSize);
    set(gca,'lineWidth',lineWidthFrame);
    set(gca,'FontName',fontName);
    set(gca,'XMinorTick','on','YMinorTick','on','ZMinorTick','on');
    set(gca,'Layer','bottom'); % Bottom/Top bottom=grid below data
    set(gca,'TickDir','out');
    
    drawnow; % required to get this hidden props
    xgrid = get(gca,'XGridHandle');
    xgrid.LineWidth = 0.5;  % default=0.5
    xgrid.GridLineStyle = '-';  % default='-'
    ygrid = get(gca,'YGridHandle');
    ygrid.LineWidth = 0.5;  % default=0.5
    ygrid.GridLineStyle = '-';  % default='-'
    zgrid = get(gca,'ZGridHandle');
    zgrid.LineWidth = 0.5;  % default=0.5
    zgrid.GridLineStyle = '-';  % default='-'
    

    %% Output
    if nargout == 0
        h = [];
    end
    
end