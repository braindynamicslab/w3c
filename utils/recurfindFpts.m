function [fixpts, Jeigs, flags] = recurfindFpts(f,initialfpts, varargin)
%RECURFINDFPTS recursively find fixed points of a dynamical system. 
%   [fixpts, Jeigs, flags] = recurfindFpts(f,initialfpts, ...)
% input:
%   f: function handle, whose zeros we like to find. 
%   initialfpts: seed fixed points, recommend using the stable steady
%   states obtained from simulation, list of fixed points, N-by-d matrix,
%   where N is number of fixed points and d the dimension of the state
%   space.
% output:
%   fixpts: list of fixed points, a N-by-d matrix.
%   Jeigs: eigenvalues of the Jacobian matrix at each fixed point.
%   flags: fsolve exit flag for each fixed point. 
% parameters:
%   fsolveOpt: options for fsolve, default optimset('fsolve')
%   maxIter: after how many levels recursion, give up searching for new
%   fixed points between a fixed list of fixed points.
%   maxFpts: stop if the maximal number of fixed points has been reached.
%   tolrFpts: the tolerance at which consider two fixed points the same,
%   i.e. they are same if abs of point-wise difference is all less than
%   this value.
%{
~ Author: Mengsen Zhang <mengsenzhang@gmail.com> 7-16-2019 ~
%}



p = inputParser;
p.addParameter('fsolveOpt',optimset('fsolve')) % parameters for fsolve
p.addParameter('maxIter',8) % maximal number of iteration when guessing the ICs
p.addParameter('maxFpts',20) % maximal number of fixed points found before stopping
p.addParameter('tolrFpts',1e-4) % threshold for considering two fixed points the same

p.parse(varargin{:})
par = p.Results;

[N_initFpts,d] = size(initialfpts);
% -- check initial fixed points
[fixpts,Jeigs,flags,idxICs]=fixedpts(f,initialfpts,par.tolrFpts,par.fsolveOpt);
foundfpts = size(uniqueTolr([initialfpts;fixpts],par.tolrFpts),1) == N_initFpts;
foundstab = all(real(Jeigs(:))<0);

if ~foundfpts
    warning('Some initial fixed points given is not a fixed point according to fsolve.')
end
if ~foundstab
    warning('Not all initial fixed points are stable.')
end

disp(['Finding all fixed points based on ' num2str(N_initFpts) ' initial fixed points.'])
[fixpts, Jeigs, flags] = findAllFpts(f,fixpts,Jeigs,flags, par);

end

function [fixpts, Jeigs, flags] = findAllFpts(f,fixpts,Jeigs,flags, par)
% find all fixed points recursively
    Nfpts = size(fixpts,1);
    fptsIdx = (1:Nfpts)';
    [fixpts_,Jeigs_,flags_,gIdx_] = findnewFpts(f,fixpts,[],1,par);
    % -- if input list of fixed points is shorter than maximum, look for
    % new ones. 
    if Nfpts<par.maxFpts && ~isempty(fixpts_)
        % -- create a new list of fixed points
        [~,fixpts_] = mixmap(fptsIdx,fixpts,gIdx_,fixpts_);
        [~,Jeigs_] = mixmap(fptsIdx,Jeigs,gIdx_,Jeigs_);
        [~,flags_] = mixmap(fptsIdx,flags,gIdx_,flags_);
        [~,uidx] = uniqueTolr(fixpts_,par.tolrFpts);
        disp(['Update list of fixed points to ' num2str(length(uidx)) ' unique points.'])
        % -- stop if list is not growing
        if length(uidx) ~= size(fixpts,1)
            fixpts = fixpts_(uidx,:);
            Jeigs = Jeigs_(uidx,:);
            flags = flags_(uidx,:);
            [fixpts, Jeigs, flags] = findAllFpts(f,fixpts,Jeigs,flags,par);
        end
    end
end

function [fixpts,Jeigs,flag,gIdx] = findnewFpts(f,oldguesses,oldgIdx,depth, par)
% search for fixed points recursively between known fixed points.
    if isempty(oldgIdx)
        oldgIdx = (1:size(oldguesses,1))';
    end
    [newguesses,newgIdx] = guessICs(oldguesses,oldgIdx);
    [fixpts,Jeigs,flag,idxICs]=fixedpts(f,newguesses,par.tolrFpts,par.fsolveOpt);
    % check if there's any new fixed points
    Nu = size(uniqueTolr([oldguesses;fixpts],par.tolrFpts),1);
    foundnew = Nu > size(oldguesses,1);
    
    if (isempty(fixpts) || (~foundnew)) && depth<par.maxIter
        % -- if no fpts found, combine old guesses, and guess again
        % inbetween
        [oldgIdx, oldguesses] = mixmap(oldgIdx,oldguesses,newgIdx,newguesses);
        [fixpts,Jeigs,flag,gIdx] = ...
            findnewFpts(f,oldguesses,oldgIdx,depth+1,par);
    else
        gIdx = newgIdx(idxICs);
    end
end

function [guesses,gIdx] = guessICs(oldguesses,oldgIdx)
    guesses = bwValues(oldguesses);
    gIdx = bwValues(oldgIdx);
end