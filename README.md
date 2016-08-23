# GRACE-filter
Contains software for filtering (destriping) GRACE Stokes coefficients

## DDK Filtering of GRACE Stokes coefficients
updated 15 October 2012 (Added DDK4 and DDK5 to package)
updated 29 January 2015 (Added DDK6 to DDK8 to package)

DDK5d9:  filtered with inverse signal degree power law 5e9*deg^4  (DDK8) weakest smoothing
DDK1d10: filtered with inverse signal degree power law 1e10*deg^4 (DDK7)         |
DDK5d10: filtered with inverse signal degree power law 5e10*deg^4 (DDK6)         |
DDK1d11: filtered with inverse signal degree power law 1e11*deg^4 (DDK5)         |
DDK5d11: filtered with inverse signal degree power law 5e11*deg^4 (DDK4)         |
DDK1d12: filtered with inverse signal degree power law 1e12*deg^4 (DDK3)         |
DDK1d13: filtered with inverse signal degree power law 1e13*deg^4 (DDK2)         |
DDK1d14: filtered with inverse signal degree power law 1e14*deg^4 (DDK1) strongest smoothing


FILTER COEFFICIENTS (in directory filtercoef):

The filter coefficients are based on a block diagonal approximation of the normal system for the month august in 2003. Together with a signal variance model which behaves as a degree dependent power law one can construct a filter matrix W which is also block diagonal:



    | W_Cl0   0        0       ..   ..   0         0      |
    |                                                     |
    | 0       WCl1     0       ..   ..   0         0      |
    |                                                     |
    | 0       0        WSl1    ..   ..   0         0      |
    |                                                     |
W = | :       :        :       \    ..   0         0      |                                                                                                         
    |                                                     |                                                                                                         
    | :       :        :       :    \    0         0      |                                                                                                         
    |                                                     |                                                                                                         
    | 0       0        0       0    0    WCl120    0      |                                                                                                         
    |                                                     |                                                                                                         
    | 0       0        0       0    0    0         WSl120 |                                                                                                         
    ------------------------------------------------------                                                                                                          
     l=2-120  l=2-120  l=2-120           l=120-120 l=120-120                                                                                                        
     ord=0    ord=1    ord=1             ord=120   ord=120                                                                                                          
     COS      COS      SIN               COS       SIN                                                                                                              
block 1       2        3                 240       241                                                                                                              
                                                                                                                                                                    
                                                                                                                                                                    
                                                                                                                                                                    
The sides of W are arranged such that orders and sine or cosine coefficients are grouped in square block diagonal matrices.                                         
                                                                                                                                                                                             
W_Clm: means an order block with Cosine coefficients of order m which varies over the degree l. Its side has size lmax-max(lmin,m)+1.                                                        
                                                                                                                                                                                             
For our case lmax=120 and lmin=2 meaning that blocks WSl120 and WCl120 have only one entry.                                                                                                  

A filtered set of a spehrical harmonic set SH can be obtained by matrix multiplication:

SH_filt= W * SH

Where the ordering scheme of the input vector SH follows that of the matrix.


## File format of DDK files
STORAGE SCHEMES and FILES:
Only the block diagonal part of the filter matrix is stored. This causes a dramatic reduction of storage space. When reading in the data one should NOT expand the matrix to its full size as this will need more than 3Gb of storage. One should therefore apply the filter blockwise.

BINARY form:

Binary Data files may be read in matlab/octave with the script read_BIN.m

On the matlab prompt type:

>> dat=read_BIN('Wbd_2-120.a_1d12p_4')


which yields:

dat =

     version: 'BINV2.1 '
        type: 'BDFULLV0'
       descr: 'ANISOTROPIC FILTER matrix with power law regularization (alpha*l^pow)           '
       nval1: 14637
       nval2: 0
       pval1: 1180123
       pval2: 0
     nblocks: 241
      ints_d: [6x24 char]
        ints: [6x1 double]
      dbls_d: [2x24 char]
        dbls: [2x1 double]
     side1_d: [14637x24 char]
    blockind: [241x1 double]
       pack1: [1180123x1 double]



The variable dat is a structure array with the following fields:
version: version of binary file type ( not relevant)
type: Type of the matrix, 'BDFULLV0' stands for a block diagonal matrix with no associated vectors

descr: Short description of the data content

nval1: size of the FULL matrix

nval2: unused yet

pval1: size of the stored matrix (packed matrix is stored in a vector with length pval1)

pval2: unused yet

nblocks: amount of diagonal blocks

ints_d: description of the integer meta data 

ints: integer meta data

dbls_d: description  of the double meta data

dbls: double meta data

side1_d: description of the elements of the side of the FULL matrix
	In case of spherical harmonics the description tag is
	TSN DEGORD or TCN DEGORD where DEG is the three digit degree and ORD is the 3 digit order of the coefficient. 
	One can retrieve the degree and order from the character array by the matlab command deg=str2num(dat.side1_d(:,4:7)) and
	ord=str2num(dat.side1_d(:,8:10))
 
	TCN denotes a cosine coefficient and TSN a sine coefficient		
	
blockind: array with indices of the last row(or column) of each block. For example:
	block n  is the submatrix: FULL(blockind(n-1)+1:blockind(n),blockind(n-1)+1:blockind(n))
	NOTE: blockind(0) is not present but should be considered zero. The blocks are stored in the order as shown in the matrix above.
	Thus block 1 corresponds to WCl0 and block 241 corresponds to WSl120

pack1: Block diagonal matrix in packed form. Blocks are stored sequentially in column major order
	so block n starts at:
	             n-1
	nstart = 1+ SUM    sz(i)^2
	             i=1

	Where the sz(i) is the size of the ith block side: 
	
	sz(i)=blockind(i)-blockind(i-1)
	

	and ends at:

	nend= nstart+ sz(n))^2	

	In matlab one can reshape the section of pack1 into a square matrix by the command:
	blockn=reshape(dat.pack1(nstart:nend),sz(n),sz(n))

As a check for a correct read one can display the  first 5 values of the pack1 matrix by issuing the following commands on the matlab prompt:

>> format long
>> dat=read_BIN('Wbd_2-120.a_1d12p_4');
>> dat.pack1(1:5)

ans =

   0.999995413476564
   0.000000054102768
   0.000000668145433
  -0.000000013374350
   0.000000077291901





REFERENCES and further reading:

For information on the theory refer to the article:
(http://www.springerlink.com/content/m75583613054m62g/)

Decorrelated GRACE time-variable gravity solutions by GFZ, and their validation using a hydrological model
Kusche, J.; Schmidt, R.; Petrovic, S.; Rietbroek, R. Journal of Geodesy, DOI: 10.1007/s00190-009-0308-3

For information on the GAC/GAA/GSM/GAD/GAB products and their use refer to the GRACE documentation:
http://isdc.gfz-potsdam.de/index.php?name=UpDownload&req=viewdownload&cid=4


