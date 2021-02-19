function [percent,tempName] = parfor_progress_mod(N,tempName,modPrint)
%PARFOR_PROGRESS Progress monitor (progress bar) that works with parfor.
%   PARFOR_PROGRESS works by creating a file called parfor_progress.txt in
%   your working directory, and then keeping track of the parfor loop's
%   progress within that file. This workaround is necessary because parfor
%   workers cannot communicate with one another so there is no simple way
%   to know which iterations have finished and which haven't.
%
%   PARFOR_PROGRESS(N) initializes the progress monitor for a set of N
%   upcoming calculations.
%
%   PARFOR_PROGRESS updates the progress inside your parfor loop and
%   displays an updated progress bar.
%
%   PARFOR_PROGRESS(0) deletes parfor_progress.txt and finalizes progress
%   bar.
%
%   To suppress output from any of these functions, just ask for a return
%   variable from the function calls, like PERCENT = PARFOR_PROGRESS which
%   returns the percentage of completion.
%
%   Example:
%
%      N = 100;
%      parfor_progress(N);
%      parfor i=1:N
%         pause(rand); % Replace with real code
%         parfor_progress;
%      end
%      parfor_progress(0);
%
%   See also PARFOR.

% By Jeremy Scheff - jdscheff@gmail.com - http://www.jeremyscheff.com/

% error(nargchk(0, 1, nargin, 'struct'));

if nargin < 1
    N = -1;
end
if ~exist('N','var') || isempty(N)
    N=1;
end
if ~exist('tempName','var') || isempty(tempName)
    [zdir,zfile] = fileparts(tempname) ;
    tempName = [ zdir '/parfor_prog_' zfile '.txt' ];
end

if nargin < 3 || isempty(modPrint)
    modPrint = 1;
end
    


percent = 0;
w = 50; % Width of progress bar

if N > 0
    f = fopen(tempName, 'w');
    if f<0
        error('Do you have write permissions for %s?', pwd);
    end
    fprintf(f, '%d\n', N); % Save N at the top of progress.txt
    fprintf(f, '%d\n', modPrint); % Save N at the top of progress.txt
    fclose(f);
    
    if nargout == 0
        disp(['  0%[>', repmat(' ', 1, w), ']']);
    end
elseif N == 0
    delete(tempName);
    percent = 100;
    
    if nargout == 0
        % disp([repmat(char(8), 1, (w+9)), char(10), '100%[', repmat('=', 1, w+1), ']']);
        fprintf(1,'Done with counter\n',percent);
    end
else
    if ~exist(tempName, 'file')
        error('parfor_progress.txt not found. Run PARFOR_PROGRESS(N) before PARFOR_PROGRESS to initialize parfor_progress.txt.');
    end
    
    f = fopen(tempName, 'a');
    fprintf(f, '1\n');
    fclose(f);
    
    f = fopen(tempName, 'r');
    progress = fscanf(f, '%d');
    fclose(f);

    lprog = length(progress)-2;    
    percent = lprog/progress(1)*100;
    modPrint = progress(2);
    
    if nargout == 0
        if modPrint > 1 && mod(lprog,modPrint) ~= 0 
           return;
        end
        perc = sprintf('%3.0f%%', percent); % 4 characters wide, percentage
        % disp([repmat(char(8), 1, (w+9)), char(10), perc, '[', repmat('=', 1, round(percent*w/100)), '>', repmat(' ', 1, w - round(percent*w/100)), ']']);
        disp(perc);
    end
end