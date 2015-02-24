function timeseries_interactive_plot(xnames)
% timeseries_interactive_plot(xnames)
%
% create interactive timeseries-plot that pops up in spm_orthviews window
%
% Tor Wager, Matthew Davidson


% find the spm window, or make one from a p-image
% -----------------------------------------------
spm_handle = findobj('Tag','Graphics');
if isempty(spm_handle) || ~ishandle(spm_handle)
    cl = pmap_threshold;
    spm_handle = gcf();
end

% setup images and mapped volumes
% -----------------------------------------------
persistent images
persistent V_images
persistent data_images

if ~exist('images', 'var') || isempty(images)
    fprintf('Setting up timeseries plot: Select images... \n');
    images = spm_get(Inf,'*img','Select images for timeseries plot');
end

if ~exist('V_images', 'var') || isempty(V_images)
    fprintf('Mapping volumes to memory... \n');
    V_images = spm_vol(images);
end

if ~exist('data_images', 'var') || isempty(data_images)
    fprintf('Loading data... \n');
    data_images = spm_read_vols(V_images);
end

% setup figure to draw to
% -----------------------------------------------
tagname = 'timeseries_interactive';

fighandle = activate_fig_window(tagname);



% set(gcf,'WindowButtonUpFcn','dat = bar_interactive_btnupfcn;')
set(spm_handle,'WindowButtonUpFcn', {@ts_interactive_callback, fighandle, V_images, data_images, tagname})

figure(spm_handle)

return




function fighandle = activate_fig_window(tagname)

    fighandle = findobj('Tag', tagname);

    if(isempty(fighandle))
        fighandle = tor_fig();
        set(fighandle, 'Tag', tagname);

        mypos = [100 100 600 200];
        set(gcf,'Position',mypos)
    else
        fighandle = fighandle(1);
        figure(fighandle);
    end

    return




% -------------------------------------------------------------------
% Callback
% -------------------------------------------------------------------

function ts_interactive_callback(spm_handle, event_data, fighandle, V_images, data_images, tagname) 
% the first two params (the source handle and event-related data) are
% required by all callbacks - see documentation
% they may be ignored if preferred; the event data may be empty

% if no images selected, don't display anything.
if isempty(V_images)
    fprintf('No images selected. Displaying nothing.\n');
    return
end

% activate window
fighandle = activate_fig_window(tagname);

% get coordinate
mm_coord = spm_orthviews('Pos');
vox_coord = mm2voxel(mm_coord',V_images(1));
mm_coord = round(mm_coord');

dat = squeeze(data_images(vox_coord(1),vox_coord(2),vox_coord(3),:));

hold off;
line_handles = plot(dat, 'k');

% if ~isempty(xnames), set(gca,'XTickLabel',xnames); end

title(sprintf('[x,y,z] = %3.0f, %3.0f, %3.0f',mm_coord(1),mm_coord(2),mm_coord(3)),'FontSize',16);


return