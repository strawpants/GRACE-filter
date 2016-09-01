%% Read in spherical harmonic coefficients from a text file
%% Usage:
%%  [out]=readSH(filename)
%%  filename (input, string): the name of the file to be opened
%%  out (output,struct): A structure containing some meta info (time, maximum degree,etc) and 
%%      out.cnm,out.snm (matrix(nmax+1,nmax+1)): Spherical harmonic coefficients in matrix form: row index and column index correspond to degree+1 and order+1 respectively
%%      out.sigcnm,out.sigsnm: Contain the errors of the coefficients with the same storage scheme (if present in the file)
%% 
%% NOTES:
%%   This routines reads in coefficients in a custom format, support for other more common formats is planned
%%
%% Copyright Roelof Rietbroek 2016
%% This software is licensed under the MIT license see https://github.com/strawpants/GRACE-filter/blob/master/LICENSE
%% URL: https://github.com/strawpants/GRACE-filter



function [out]=readSH(filename)

fid=fopen(filename,'r');

%% get the first line with meta info

meta=textscan(fgetl(fid),'%s %d %f %f %f %f\n');
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



%put the goodies in matrix
out.cnm=zeros(out.nmax+1,out.nmax+1);
out.snm=zeros(out.nmax+1,out.nmax+1);
out.cnm(sub2ind([out.nmax+1,out.nmax+1],dat{1}+1,dat{2}+1))=dat{3};
out.snm(sub2ind([out.nmax+1,out.nmax+1],dat{1}+1,dat{2}+1))=dat{4};

%same for the errors
out.sigcnm=zeros(out.nmax+1,out.nmax+1);
out.sigsnm=zeros(out.nmax+1,out.nmax+1);

if(sigmas)
  out.sigcnm(sub2ind([out.nmax+1,out.nmax+1],dat{1}+1,dat{2}+1))=dat{5};
  out.sigsnm(sub2ind([out.nmax+1,out.nmax+1],dat{1}+1,dat{2}+1))=dat{6};
end


fclose(fid);



end

