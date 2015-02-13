% This program calculates the aberration dependence of contrast in DIC microscopy
% Copyright (C) 2015 Moritz Esslinger
% moritz.esslinger@web.de
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
% GNU General Public License for more details.
% You should have received a copy of the GNU General Public License along
% with this program; if not, write to the Free Software Foundation, Inc.,
% 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
% 



function C = matrix_multiplication( A, B, summation_dimA, summation_dimB)
% Multiplies Arrays like a matrix multiplication. 
% Dimensionality may be higher than 2D.
%
% Input: A              : Array A. At least one-dimensional.
%        B              : Array B. At least one-dimensional.
%        summation_dimA : dimension (index number) over which 
%                         the summation is performed for array A.
%        summation_dimB : dimension (index number) over which 
%                         the summation is performed for array B.
%
%
% Output: C             : result of the matrix multiplication C = A * B.
%                         For a N-dimensional matrix A and a
%                         M-dimensional matrix B, the result C is
%                         (N - 1) + (M - 1) dimensional. Ordering of
%                         indices of C is first all indices of A and then all
%                         indices of B, except for those over which the
%                         summation took place.
%
% Examples:
%           Scalar product of two vectors
%             C = [1,2,3] * [4,5,6]'
%             C = matrix_multiplication( [1,2,3] , [4,5,6] , 1 , 1 )
%
%           2D times 2D Matrix Multiplication
%             A       = zeros(4,5);
%             A(1:20) = 1:20;
%             B       = zeros(5,7);
%             B(1:35) = 42 + 1:35;
%             C = A * B
%             C = matrix_multiplication(A,B,2,1)
%
%           3D times 3D matrix multiplication
%             A_ijk  = zeros(2,3,4);
%             B_kmn  = zeros(4,7,5);
%             C_ijmn = matrix_multiplication(A_ijk, B_kmn, 3, 1)
%
%           3D times 3D matrix multiplication
%             A_ijk  = zeros(2,3,4);
%             B_ljn  = zeros(4,7,5);
%             C_ikln = matrix_multiplication(A_ijk, B_ljn, 2, 2)
%
%


sa = size(A);
summation_lengtha = prod(sa(summation_dimA));
nonsummation_dimA = 1:length(sa);
nonsummation_dimA(summation_dimA) = [];
nonsummation_lengtha = numel(A) / summation_lengtha;

sb = size(B);
summation_lengthb = prod(sb(summation_dimB));
nonsummation_dimB = 1:length(sb);
nonsummation_dimB(summation_dimB) = [];
nonsummation_lengthb = numel(B) / summation_lengthb;

if summation_lengtha ~= summation_lengthb
    disp('ERROR: Inner matrix dimensions must agree')
    disp([ '    size(A) = [ ', num2str(sa), ' ] , summation over axis ', num2str(summation_dimA) ])
    disp([ '    size(B) = [ ', num2str(sb), ' ] , summation over axis ', num2str(summation_dimB) ])
    disp([ '    ', num2str(summation_lengtha), ' ~= ', num2str( summation_lengthb)])
else
    
    % reshape the input into two-dimensional matrices, 
    Atilde = reshape( permute(A, [nonsummation_dimA, summation_dimA]),    nonsummation_lengtha, summation_lengtha);    
    Btilde = reshape( permute(B, [summation_dimB, nonsummation_dimB]),    summation_lengthb, nonsummation_lengthb);

    % the actual multiplication
    Ctilde = Atilde * Btilde;
    
    % reshape back
    % First  index Ctilde(i, ... ) resembles indices nonsummation_dimA of A
    % Second index Ctilde(... , i) resembles indices nonsummation_dimB of B
    
    % size of the final array C
    sc = [ sa(nonsummation_dimA), sb(nonsummation_dimB) ];
    C = reshape( Ctilde, sc);
    
end
