function [h] = nicearea(varargin)
    
    
    % Parameters
    fontSize       = 25;
    lineWidthFrame = 1;
    lineWidth      = 1;
    fontName       = 'Times';
    
    
    %% If there is no specification of lineWidth, use default of lineWidth
    lineWidth_flag = 0;
    for i=1:1:nargin
        if strcmpi('lineWidth',varargin{i})
            lineWidth_flag = 1;
        end
    end
    if ~lineWidth_flag
        varargin{nargin+1}='lineWidth';
        varargin{nargin+2}=lineWidth;
    end
    
    
    %% Plot
    h=area(varargin{:});
    
    
    %% Make nicer
    box on;
    grid on;
    set(gca,'fontSize',fontSize);
    set(gca,'lineWidth',lineWidthFrame);
    set(gca,'FontName',fontName);
    set(gca,'XMinorTick','on','YMinorTick','on');
    ax = gca; ax.Layer = 'top'; % Grid in front of data
    
    %% Output
    if nargout == 0
        h = [];
    end
end