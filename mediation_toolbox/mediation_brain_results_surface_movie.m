function mov = mediation_brain_results_surface_movie
%
% Make a movie of mediation_brain_results surfaces
% mov = mediation_brain_results_surface_movie;
%
% Creates mymovie.avi in current directory
%
% Steps:
% 1) get clusters from mediation_brain_results
% 2) make surface figures
% 3) run this function
%
% Example code:
%
% [clpos, clneg, clpos_data, clneg_data] = mediation_brain_results('rob0', 'thresh', [.000001 .00001 .0008], 'size', [1 1 1], ...
% mediation_brain_surface_figs(clpos, clneg);
% mov = mediation_brain_results_surface_movie;
%
% What can go wrong?
% Among other things, errors with OpenGL rendering in matlab,
% and/or movie saving errors.  Not sure why this happens,
% but you may want to make sure you can save .png images
% or other images of the figure before you try making the 
% movie.  
% With some versions of matlab, longer movies have not worked
% well for me, but I haven't figured out why.
%
% may want to use getframe matlab movie style instead
% change movie_tools code to FRAMESTYLE = 'matlab' then run
%
% Tor Wager, June 2008

% Get handles
% -------------------
f1 = create_figure('Movie_fig');
f1h = gca;

forig = create_figure('Mediation Surface Figures', 1, 1, 1);
surfhan = findobj(forig, 'Type', 'Patch');

% Copy surface
% -------------------
sh_new = copyobj(surfhan(1), f1h);

figure(f1);
axes(f1h);

lighting gouraud; axis off; axis image; lightRestoreSingle;
drawnow

% Make movie
% -------------------
mov = movie_tools('batch360.1');

sh_new_subctx = copyobj(surfhan([2 4:end-3]), f1h);

mov = movie_tools('transparent', 1, 0, sh_new, mov, 2.5);
mov = movie_tools('rotate', 90, 0, mov);

mytype = class(mov);
if strcmp(mytype, 'avifile')
    % we already closed the movie in movie_tools: avifile
else
    % we need to write the AVI file
    movie2avi(mov, 'surface_360.avi', 'fps', 5);
end

end

