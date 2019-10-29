%% Test routine to see if the anisotropic filter routine for GRACE spherical harmonic coefficients is implemented
%% correctly in matlab and octave
%% The script contains three experiments:
%% Experiment 1: filter a set of spherical harmonic coefficients which have a maximum degree which equals that of the filter
%% Experiment 2: filter a set of spherical harmonic coefficients which have a maximum degree larger than that of the filter
%%               The desired behavior is that the filtered coefficients will have zero values where their degrees exceedt that of the matrix
%% Experiment 3: filter a set of spherical harmonic coefficients which have a maximum degree which is smaller than that of the filter
%% Experiment 4: Check the error propagation


%% NOTES: 
%%       One should always FILTER THE STOKES RESIDUALS wrt. a static gravity field. The provided test coefficients are already residuals
%%       please consult the README.md in the root of the git repo for more information
%%
%%       In the script below, the real spherical harmonic coefficients are gathered in (lower triangular) matrices where rows correspond to the degree-1 and the column to the order-1
%%          Examples:
%%              cnm(5,2) and snm(2,1) corresponds to C41 and S10 respectively
%%              cNM=cnm(N+1,M+1)
%%              cnm(2,6)=0 as the order exceeds the degree

%% Copyright Roelof Rietbroek 2016
%% This software is licensed under the MIT license see https://github.com/strawpants/GRACE-filter/blob/master/LICENSE
%% URL: https://github.com/strawpants/GRACE-filter


clear all;

%% Here's where we find some additional matlab scripts needed
addpath('../src/matlab');

%% Read the block diagonal filter matrix (the file below corresponds to DDK2 and has filter coefficients from degrees 2 to 120)
Wbd=read_BIN('../data/DDK/Wbd_2-120.a_1d13p_4');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Experiment 1 (nmaxinput== nmaxfilter)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp("Experiment 1: testing up to maximum filter degree");

%% read in input coefficients of maximum degree and order 120
input=readSH('GSM-2_2008122-2008153_0030_EIGEN_G---_0004in');


%% Filter the input coefficients
[cnmfilt,snmfilt]=filterSH(Wbd,input.cnm,input.snm);

%read in independently filtered values for validation
checkvals=readSH('GSM-2_2008122-2008153_0030_EIGEN_G---_0004out');

%check if the differences are neglible (numerically insignificant)
if (max(max(abs(checkvals.cnm-cnmfilt))) > eps('double'))
   error('Experiment 1: filtered cnm matrix not the same as independent values'); 
end

if (max(max(abs(checkvals.snm-snmfilt))) > eps('double'))
   error('Experiment 1: filtered snm matrix not the same not the same as independent values'); 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Experiment 2 (nmaxinput > nmaxfilter)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp("Experiment 2: testing with a dataset which has larger maximum degree than the filter");
%% extend the data by this amount of degrees:
nmaxadd=5;

%% Extend the input coefficients of the first experiment set coefficients to one where n > nmaxfilter
cnm=ones(input.nmax+1+nmaxadd,input.nmax+1+nmaxadd);
snm=ones(input.nmax+1+nmaxadd,input.nmax+1+nmaxadd);
cnm(1:input.nmax+1,1:input.nmax+1)=input.cnm;
snm(1:input.nmax+1,1:input.nmax+1)=input.snm;

[cnmfilt,snmfilt]=filterSH(Wbd,cnm,snm);

%% Also extend the validation coefficients and set coefficients to zero where n > nmaxfilter (this is the desired behavior of the filter routine)
checkcnm=zeros(input.nmax+1+nmaxadd,input.nmax+1+nmaxadd);
checksnm=zeros(input.nmax+1+nmaxadd,input.nmax+1+nmaxadd);
checkcnm(1:input.nmax+1,1:input.nmax+1)=checkvals.cnm;
checksnm(1:input.nmax+1,1:input.nmax+1)=checkvals.snm;

%% validate the filtered coefficients
if (max(max(abs(checkcnm-cnmfilt))) > eps('double'))
   error('Experiment 2: filtered cnm matrix not the same as independent values'); 
end

if (max(max(abs(checksnm-snmfilt))) > eps('double'))
     error('Experiment 2: filtered snm matrix not the same as independent values'); 
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Experiment 3 (nmaxinput < nmaxfilter)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp("Experiment 3: testing with a dataset which has a smaller maximum degree than the filter");
%% Read in the truncated (nmax=60) filter coefficients
input=[];
input=readSH('GSM-2_2008122-2008153_0030_EIGEN_G---_0004lmax60in');
%% filter coefficients
[cnmfilt,snmfilt]=filterSH(Wbd,input.cnm,input.snm);

%% read in indpendent filter results of the truncated version
checkvals=readSH('GSM-2_2008122-2008153_0030_EIGEN_G---_0004lmax60out');

if (max(max(abs(checkvals.cnm-cnmfilt))) > eps('double'))
   error('Experiment 3: filtered cnm matrix not the same as independent values'); 
   
end

if (max(max(abs(checkvals.snm-snmfilt))) > eps('double'))
   error('Experiment 3: filtered snm matrix not the same as independent values'); 
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Experiment 4  Test error propagation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp("Experiment 4: Testing error propagation");
%% read in input coefficients of maximum degree and order 120
input=readSH('GSM-2_2008122-2008153_0030_EIGEN_G---_0004in');


%% Filter the input coefficients and propagate standard error
[cnmfilt,snmfilt,sigcnmfilt,sigsnmfilt]=filterSH(Wbd,input.cnm,input.snm,input.sigcnm,input.sigsnm);

%read in independently filtered values for validation
checkvals=readSH('GSM-2_2008122-2008153_0030_EIGEN_G---_0004out');

%check if the differences are neglible (numerically insignificant)
if (max(max(abs(checkvals.sigcnm-sigcnmfilt))) > eps('double'))
   error('Experiment 1: filtered sigcnm matrix not the same as independent values'); 
end

if (max(max(abs(checkvals.sigsnm-sigsnmfilt))) > eps('double'))
   error('Experiment 1: filtered sigsnm matrix not the same not the same as independent values'); 
end
%% If we end up here we did all tests succesfully

fprintf('Filter routines successfully completed\n');
