%% Read in spherical harmonic coefficients from a text file
%% Usage:
%%  [out]=readSH(filename,varargin)
%%
%%  filename (input, string): the name of the file to be opened
%%  varargin can be pairs of the following:
%%      'form',X : set output format of the coefficients
%%      X=0: Matrix format (default)
%%      X=1: Vectorized format 
%%
%%      'nmax',NMAX : restrict output to a lesser degree than file
%% 
%%
%%  output when form=0 (matrix): 
%%  out (output,struct): A structure containing some meta info (time, maximum degree,etc) and 
%%      out.cnm,out.snm (matrix(nmax+1,nmax+1)): Spherical harmonic coefficients in matrix form: row index and column index correspond to degree+1 and order+1 respectively
%%      out.sigcnm,out.sigsnm: Contain the errors of the coefficients with the same storage scheme (if present in the file)
%% 
%%  output when form=1 (vector)
%%  out is a structure as above but
%%  out.cnm is a vector nmax*(nmax+1)
%%  out.n is the corresponding degree
%%  out.m is the corresponding order
%%  out.descr is a 24 character description
%% NOTES:
%%   This routines reads in coefficients in a custom format, support for other more common formats is planned
%%
%% Copyright Roelof Rietbroek 2017
%% This software is licensed under the MIT license see https://github.com/strawpants/GRACE-filter/blob/master/LICENSE
%% URL: https://github.com/strawpants/GRACE-filter



function [out]=readSH(filename,varargin)

format=0; %0: 'matrix', 1: 'vector'
nmax=-1;

i=1;
while i<nargin 
    
    switch varargin{i}
    case 'form'
        format=varargin{i+1};
        i=i+1;
    case 'nmax'
        nmax=varargin{i+1};
        i=i+1;
    end
    i=i+1;
end

if format > 1
    error('output form not understood');
end 


fid=fopen(filename,'r');

%% get the first line with meta info
meta=textscan(fgetl(fid),'%s %d %f %f %f %f\n');
%meta=textscan(fgetl(fid),'%s %d %f %f %f %f\n');
out.nmax=meta{2};
out.yearstart=meta{3};
out.yearcenter=meta{4};
out.yearend=meta{5};



%try to find out how many columns are present by scanning the first data line
fposition=ftell(fid);
testl=textscan(fgetl(fid),'%d %d %f %f %f %f\n');
if(isempty(testl{5}))
  sigmas=0;
else
  sigmas=1;
end

%% rewind the file to the beginning of the data (after the first meta info line)
fseek(fid,fposition,'bof');

if (sigmas)
  %read the remainder in a temporary matrix with sigmas
  dat=textscan(fid,'%d %d %f %f %f %f\n');
 else
  %read the remainder in a temporary matrix without sigmas
  dat=textscan(fid,'%d %d %f %f\n');
end

fclose(fid);


%%Check if output degree is going to be correctly restricted
if nmax>=0
    if nmax > out.nmax
        error(['File supports degrees up to ' num2str(out.nmax) ', requested:' num2str(nmax)])
    end
    
    out.nmax=nmax;

    %restrict input
    indxv=find(le(dat{1},nmax));
    dat{1}=dat{1}(indxv);
    dat{2}=dat{2}(indxv);
    dat{3}=dat{3}(indxv);
    dat{4}=dat{4}(indxv);
    if sigmas
        dat{5}=dat{5}(indxv);
        dat{6}=dat{6}(indxv);

    end
else
    nmax=out.nmax;
end 

switch format
case 0 %matrix
    %put the goodies in matrix
    out.cnm=zeros(nmax+1,nmax+1);
    out.snm=zeros(nmax+1,nmax+1);
    out.cnm(sub2ind([nmax+1,nmax+1],dat{1}+1,dat{2}+1))=dat{3};
    out.snm(sub2ind([nmax+1,nmax+1],dat{1}+1,dat{2}+1))=dat{4};

    %same for the errors
    out.sigcnm=zeros(nmax+1,nmax+1);
    out.sigsnm=zeros(nmax+1,nmax+1);

    if(sigmas)
      out.sigcnm(sub2ind([nmax+1,nmax+1],dat{1}+1,dat{2}+1))=dat{5};
      out.sigsnm(sub2ind([nmax+1,nmax+1],dat{1}+1,dat{2}+1))=dat{6};
    end
case 1 %vectorized format 

    cdescr=sprintf('GCN %3d%3d              ',[dat{1} dat{2}]');
    out.cnm=dat{3};

    if sigmas
        out.sigcnm=dat{5};
    end
    out.nm=[dat{1} dat{2}]; 
    %squeeze out m=0 sine coefficients
    indxv=find(ne(dat{2},0));
    dat{1}=dat{1}(indxv);
    dat{2}=dat{2}(indxv);
    dat{4}=dat{4}(indxv);
    sdescr=sprintf('GSN %3d%3d              ',[dat{1} dat{2}]');
    out.nm=[out.nm; dat{1} dat{2}];
    out.cnm=[out.cnm; dat{4}];


    if sigmas
        out.sigcnm=[out.sigcnm; dat{6}];
    end

    out.descr=reshape([cdescr sdescr],24,size(out.cnm,1))';



end



end

