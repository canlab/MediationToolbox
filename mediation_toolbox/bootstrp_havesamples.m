function [bootstat, bootsam] = bootstrp_havesamples(bootsam,bootfun,varargin)
%function [bootstat, bootsam] = bootstrp_havesamples(bootsam,bootfun,varargin)
%
% Tor Wager, Sept. 2009
%
% Drop-in replacement for Matlab's bootstrp, using pre-specified bootstrap samples
%
% Initialize matrix to identify scalar arguments to bootfun.
la = length(varargin);
scalard = zeros(la,1);

% find out the size information in varargin.
n = 1;
for k = 1:la
   [row,col] = size(varargin{k});
   if max(row,col) == 1
      scalard(k) = 1;
   end
   if row == 1 && col ~= 1
      row = col;
      varargin{k} = varargin{k}(:);
   end
   n = max(n,row);
end

if isempty(bootfun)
   bootstat = zeros(nboot,0);
   return
end

haveallsamples = 1;
nboot = size(bootsam, 2);

% Get result of bootfun on actual data and find its size.
bootstat = feval(bootfun,varargin{:});

% Initialize an array to contain the results of all the bootstrap
% calculations, preserving the output type
bootstat(nboot,1:numel(bootstat)) = bootstat(:)';

% Do bootfun - nboot times.
if la==1 && ~haveallsamples && ~any(scalard)
   % For special case of one non-scalar argument and one output, try to be fast
   X1 = varargin{1};
   for bootiter = 1:nboot
      onesample = ceil(n*rand(n,1));
      tmp = feval(bootfun,X1(onesample,:));
      bootstat(bootiter,:) = (tmp(:))';
   end
elseif la==2 && ~haveallsamples && ~any(scalard)
   % For two non-scalar arguments and one output, try to be fast
   X1 = varargin{1};
   X2 = varargin{2};
   for bootiter = 1:nboot
      onesample = ceil(n*rand(n,1));
      tmp = feval(bootfun,X1(onesample,:),X2(onesample,:));
      bootstat(bootiter,:) = (tmp(:))';
   end
else
   % General case
   db = cell(la,1);
   for bootiter = 1:nboot
      if haveallsamples
         onesample = bootsam(:,bootiter);
      else
         onesample = ceil(n*rand(n,1));
      end
      for k = 1:la
         if scalard(k) == 0
            db{k} = varargin{k}(onesample,:);
         else
            db{k} = varargin{k};
         end
      end
      tmp = feval(bootfun,db{:});
      bootstat(bootiter,:) = (tmp(:))';
   end
end