%% Function which reads in a binary file containing symmetric/full or block
%% diagonal matrices and associated vectors and parameters
%% 
%% Usage: dat=read_BIN(file)
%%      file (input,string): input file
%%      dat (output,struct): Structure with the file content
%%      the matrix remains in packed form (dat.pack1 field)
%%
%% Or: dat=read_BIN(file,'F')
%%      also expands the matrix to its full form (dat.mat1 field)
%%      Warning: this option may cause excessive RAM memory use with large matrices
%% 
%% Initial version: 7-1-2008
%% Updated: 29-07-2008: function now also works with Octave
%% Updated 6-08-08: new Binary version can also be read in
%% Updated 5-12-08: unpacking is perfomed with a compiled mex fortran routine
%% Updated 01-08-2010 allows version 2.5 to be read
%%
%% Copyright Roelof Rietbroek 2016
%% This software is licensed under the MIT license see https://github.com/strawpants/GRACE-filter/blob/master/LICENSE
%% URL: https://github.com/strawpants/GRACE-filter
function dat=read_BIN(file,varargin)

unpack=false;


for i=1:size(varargin,2)
  switch varargin{i}
   case {'F'}
    unpack=true; % unpack matrix in full size
  
  end
end



%open file for read acces in little endian format
[fid,message]=fopen(file,'r','ieee-le');
%check for errors
 if (fid == -1) 
   message
   dat=[];
   return; 
 end

%read endian checker
endian=fread(fid,1,'uint16')';
 
if endian ~= 18754 % reopen file in big endian mode
  fclose(fid);
  [fid,message]=fopen(file,'r','ieee-be');
  endian=fread(fid,1,'uint16')';
end

%read BINARY version and type from file

dat.version(1,3:8)=fread(fid,6,'uint8=>char')';
dat.version(1,1:2)='BI';
% %convert character version to a number (needed for compatibility with
% %older versions)
dat.ver=str2num(dat.version(1,5:8));

dat.type(1,1:8)=fread(fid,8,'uint8=>char')';
dat.descr(1,1:80)=fread(fid,80,'uint8=>char')';

%read indices
%integers:inint,indbls,inval1,inval2,ipval1,ipval2
metaint=fread(fid,4,'uint32');
%put index data in structure array
dat.nints=metaint(1);
dat.ndbls=metaint(2);
dat.nval1=metaint(3);
dat.nval2=metaint(4);


if dat.ver < 2.4 % compatibility clause ( pval are long integers in newer versions)
  metaint=fread(fid,2,'uint32');
  dat.pval1=metaint(1);
  dat.pval2=metaint(2);
else
  metaint=fread(fid,2,'uint64');
  dat.pval1=metaint(1);
  dat.pval2=metaint(2);
end



if dat.ver <= 2.1 %compatibility clause
  switch dat.type
   case {'SYMV0___','BDFULLV0','BDSYMV0','BDFULLVN'}
    dat.nvec=0;
    dat.pval2=1;
   case {'SYMV1___'}
    dat.nvec=1;
    dat.pval2=1;
   case {'SYMV2___'}
    dat.nvec=2;
    dat.pval2=1;
   case{'FULLSQV0'}
    dat.nvec=0;
    dat.pval2=dat.pval1;    
  end
  dat.nread=0;  
  dat.nval2=dat.nval1;
else
  dat.nvec=fread(fid,1,'integer*4');
  dat.nread=fread(fid,1,'integer*4');  
end

%Type dependent index data
switch dat.type
 case {'BDSYMV0_','BDFULLV0','BDSYMVN_','BDFULLVN'}
  %read additional nblocks parameter
  dat.nblocks=fread(fid,1,'integer*4');
end 

%get meta data 
%readme array
if(dat.nread > 0)
  list=fread(fid,dat.nread*80,'uint8=>char');
  dat.readme=reshape(list,80,dat.nread)';
end

%integers
if(dat.nints > 0)
  if dat.ver <= 2.4  %compatibility clause
    list=fread(fid,dat.nints*24,'uint8=>char');
    dat.ints_d=reshape(list,24,dat.nints)';
    dat.ints=fread(fid,dat.nints,'integer*4');
  else
    list=fread(fid,dat.nints*24,'uint8=>char');
    dat.ints_d=reshape(list,24,dat.nints)';
     dat.ints=fread(fid,dat.nints,'int64');
     
  end
end


%doubles
if(dat.ndbls > 0)
  list=fread(fid,dat.ndbls*24,'uint8=>char');
  dat.dbls_d=reshape(list,24,dat.ndbls)';
  dat.dbls=fread(fid,dat.ndbls,'real*8');
end




%side description meta data
list=fread(fid,dat.nval1*24,'uint8=>char');
%reshape characters and put in dat struct array
dat.side1_d=reshape(list,24,dat.nval1)';



%type specific meta data
switch dat.type
 case {'BDSYMV0_','BDFULLV0','BDSYMVN_','BDFULLVN'}
  %read additional nblocks parameter
  dat.blockind=fread(fid,dat.nblocks,'integer*4');

end

switch  dat.type
 case{'BDFULLV0','BDFULLVN','FULLSQV0','FULLSQVN'}

  if dat.ver <=2.2 % compatibility clause
    dat.side2_d=dat.side1_d;
  else
    list=fread(fid,dat.nval2*24,'uint8=>char');
    %reshape characters and put in dat struct array
    dat.side2_d=reshape(list,24,dat.nval2)';
  end

 case{'FULL2DVN'}
  list=fread(fid,dat.nval2*24,'uint8=>char');
  %reshape characters and put in dat struct array
  dat.side2_d=reshape(list,24,dat.nval2)';
end 


%data (type dependent)

%vectors
if dat.nvec >0
  for i=1:dat.nvec
    dat.vec(:,i)=fread(fid,dat.nval1,'real*8');
  end
end


%read matrix data
dat.pack1=fread(fid,dat.pval1*dat.pval2,'real*8');

%CLOSE FILE
fclose(fid);


if(~unpack)
  return 
end
%unpack if requested

switch dat.type
  
 case {'SYMV0___','SYMV1___','SYMV2___','SYMVN___'}
  
  %copy data from packed vector to full array
  %slow version with MATLAB loops:
  st=0;
  for i=1:dat.nval1
    dat.mat1(1:i,i)=dat.pack1(st+1:st+i);
    st=st+i;
  end
  %%end slow version
  
  %quick version with mex function
%  dat.mat1=symunpack(dat.pack1);
  %end quick version

  dat=rmfield(dat,'pack1');

 case {'BDSYMV0_','BDSYMVN_'}
  %fill first block
  sz=dat.blockind(1);
  shift1=1;
  shift2=sz*(sz+1)/2;
  dat.mat1=symunpack(dat.pack1(shift1:shift2));
  %loop over blocks
  
  shift1=shift2;
  for i=2:dat.nblocks
    sz=dat.blockind(i)-dat.blockind(i-1);
    shift2=shift1+sz*(sz+1)/2;
    dat.mat1=blkdiag(dat.mat1,symunpack(dat.pack1(shift1+1:shift2)));
    shift1=shift2;
  end

  dat=rmfield(dat,'pack1');
 case{'BDFULLV0','BDFULLVN'}
  
  %fill first block
  sz=dat.blockind(1);
  shift1=1;
  shift2=sz^2;
  dat.mat1=reshape(dat.pack1(shift1:shift2),sz,sz);
  %loop over blocks
  
  shift1=shift2;
  for i=2:dat.nblocks
    sz=dat.blockind(i)-dat.blockind(i-1);
    shift2=shift1+sz^2;
    dat.mat1=blkdiag(dat.mat1,reshape(dat.pack1(shift1+1:shift2),sz,sz));
    shift1=shift2;
  end

  dat=rmfield(dat,'pack1');
  

 case{'FULLSQV0','FULLSQVN','FULL2DVN'}

    dat.mat1=reshape(dat.pack1,dat.nval1,dat.nval2);
      dat=rmfield(dat,'pack1');

end
  
end 


