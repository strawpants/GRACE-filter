%% Filter spherical harmonic coefficients witha block diagonal fitler matrix
%% Usage:
%%  [cnmfilt,snmfilt]=filterSH(W,cnm,snm)
%%  [cnmfilt,snmfilt,sigcnmfilt,sigsnmfilt]=filterSH(W,cnm,snm,sigcnm,sigsnm) (also propagate standard deviations)
%%  W (input, struct): Struct containing the filter matrix (read by BIN_readslow)
%%  cnm,snm (input,matrix(nmax+1,nmax+1)): Spherical harmonic coefficients in matrix form: row index and column index correspond to degree+1 and order+1 respectively
%%  cnmfilt,snmfilt (output,matrix(nmax+1,nmax+1)): Filtered spherical harmonic coefficients
%% optional arguments and output:
%% sigcnm,sigsnm (input, matrix(nmax+1,nmax+1)): standard errors of the input coefficients
%% sigcnmfilt,sigsnmfilt (output, matrix(nmax+1,nmax+1)): standard errors of the filtered output coefficients(using diagonal error propagation)
%%
%% TODO: extend the function to also allow error propagation (using varargin for input)
%% Copyright Roelof Rietbroek 2019
%% This software is licensed under the MIT license see https://github.com/strawpants/GRACE-filter/blob/master/LICENSE
%% URL: https://github.com/strawpants/GRACE-filter
function [varargout]=filterSH(varargin)
    %extract filter matrix 
    W=varargin{1};
    cnm=varargin{2};
    snm=varargin{3};
    properror=0;

    if nargin == 5
        sigcnm=varargin{4};
        sigsnm=varargin{5};
        properror=1;
    end
    if (nargin != 3) && (nargin !=5)
        error("either call the function as [cnmfilt,snmfilt]=filterSH(W,cnm,snm) or [cnmfilt,snmfilt,sigcnmfilt,sigsnmfilt]=filterSH(W,cnm,snm,sigcnm,sigsnm)");
    end



%% check if we have a block diagonal filter matrix
if(isfield(W,'type'))
    if(strcmp(W.type,'BDFULLV0'))
        %% maximum degree of the input coefficients
        nmax=size(cnm,1)-1;
        
        %% Extract the minimum and maximum degree supported by the filter matrix
        for i=1:size(W.ints_d,1)
            if(strncmp(W.ints_d(i,:),'Lmax',4))
                nmaxfilt=W.ints(i);
            end
            if(strncmp(W.ints_d(i,:),'Lmin',4))
                nminfilt=W.ints(i);
            end
        end
        
        %% Determine the output maximum degree (limited by either the filter or input data)
        nmaxout=min(nmax,nmaxfilt);
    
        %% Reserve space for output (will have same size as input) and set to zero
        cnmfilt=zeros(nmax+1,nmax+1);
        snmfilt=zeros(nmax+1,nmax+1);
    
        if properror
            sigcnmfilt=zeros(nmax+1,nmax+1);
            sigsnmfilt=zeros(nmax+1,nmax+1);
        end

        %% Loop parameter indicating the previous block number
        lastblckind=0;
        %% Loop parameter indicating the end position in the packed matrix of the previous block
        lastindex=0;
        
        %% loop over the available blocks
        for iblk=1:W.nblocks;
            %% Get the order of the block from the block index
            order=floor(iblk/2);
            
            %% Break loop if the degrees of the block are larger than the degrees of the input
            if order > nmaxout
                break;
            end
            
            
            %% get the trigonometric part type (cosine=1 sin=0) from the index
            trig=floor(mod(iblk+gt(iblk,1),2));
            
            % Compute the size of the side of the stored block
            sz=W.blockind(iblk)-lastblckind;
            
            %% Initialize the filter order block to a unit diagonal matrix
            %% But accomodate space in the block also for elements with smaller
            %% degrees than nminfilt
            blockn=diag(ones(nmaxfilt+1-order,1));
           
            %% Minimum (stored) degree for this particular block (may be limited by the mininum degree supported by the filter)
            nminblk=max(nminfilt,order);
            %% This may cause an offset which we need to account for when expanding the block in a fully occupied order blocks
            shft=nminblk-order+1;
            
            %% unpack the stored filterblock (vector) in a fully occupied orderblock matrix
            blockn(shft:end,shft:end)=reshape(W.pack1(lastindex+1:lastindex+sz^2),sz,sz);
            
            %% Filter the input coefficients (this is in fact just a matrix vector multiplication)
            if(trig)
                cnmfilt(order+1:nmaxout+1,order+1)=blockn(1:nmaxout+1-order,1:nmaxout+1-order)*cnm(order+1:nmaxout+1,order+1);
            else
                snmfilt(order+1:nmaxout+1,order+1)=blockn(1:nmaxout+1-order,1:nmaxout+1-order)*snm(order+1:nmaxout+1,order+1);
            end
            if properror 
                blockn=blockn.^2;
                if(trig)
                    sigcnmfilt(order+1:nmaxout+1,order+1)=sqrt(blockn(1:nmaxout+1-order,1:nmaxout+1-order)*(sigcnm(order+1:nmaxout+1,order+1).^2));
                else
                    sigsnmfilt(order+1:nmaxout+1,order+1)=sqrt(blockn(1:nmaxout+1-order,1:nmaxout+1-order)*(sigsnm(order+1:nmaxout+1,order+1).^2));
                end
            end
            %% Prepare the loop variables for next block
            lastblckind=W.blockind(iblk);
            lastindex=lastindex+sz^2;
        end

        varargout{1}=cnmfilt;
        varargout{2}=snmfilt;
        if properror
            varargout{3}=sigcnmfilt;
            varargout{4}=sigsnmfilt;
        end
    else
        error('Not a block diagonal matrix');
    end
else
    error('Unknown filter type');
end


end
