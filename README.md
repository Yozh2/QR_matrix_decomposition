# QR_matrix_decomposition
Using Householder method

## Main source file: ``` qr.c ```
## Purpose	
Demonstrate the QR decomposition. The method of Householder reflections should be used.

## Input
rectangular MxN matrix A from the file ```A.txt``` which is the matrix of linear algebraic equasions system. The input file contains the ```M and N``` numbers at the first line. Below them there are the coefficients of the matrix.

## Output
the `Q` matrix, the `R` matrix and the `QR` multiplication (`A` matrix)

## Compilation
Using gcc:
```
gcc -g -Wall -O3 -o qr qr.c
```
or via makefile:
```
make
```

## Usage    
```
./qr
```
The file `A.txt` containing the input matrix should be in the same folder with the executable file

### Example
Input file `A.txt` contains:

```
  5 3
  12.000  -51.000    4.000
   6.000  167.000  -68.000
  -4.000   24.000  -41.000
  -1.000    1.000   -0.000
   2.000   -0.000    3.000
```
after running `./qr` we will get output
```
    Q:
    0.846   -0.391    0.343    0.082    0.078
    0.423    0.904   -0.029    0.026    0.045
   -0.282    0.170    0.933   -0.047   -0.137
   -0.071    0.014   -0.001    0.980   -0.184
    0.141   -0.017   -0.106   -0.171   -0.969

    R:
   14.177   20.667  -13.402
   -0.000  175.043  -70.080
    0.000    0.000  -35.202
   -0.000   -0.000   -0.000
    0.000    0.000   -0.000

    Q * R:
   12.000  -51.000    4.000
    6.000  167.000  -68.000
   -4.000   24.000  -41.000
   -1.000    1.000   -0.000
    2.000   -0.000    3.000
```
## Licence
This code is free for any use and anyone can modify it as he wants.
