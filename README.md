# Flex 4.0

## Table of contents
1. [General description](#General-description)
2. [Installing suggested](#Installing-suggested)
3. [Input files](#Input-files)
4. [Output files](#Output-files) 

## General description
Implementation of the numerical lattice Self-Consistent Field (SCF) method for simulations of binary brushes in selective solvent.

A detailed description of the method and algorithm is given in the article ["Chains Stiffness Effect on the Vertical Segregation of Mixed Polymer Brushes in Selective Solvent"](https://www.mdpi.com/2073-4360/15/3/644/ "Polymers 2023, 15(3), 644").

## Installing suggested 
To install the calculation package **Flex_2_2**, download this resource, go to the **Source** folder and execute the **make** command in terminal.
The code will be assembled and compiled according to the instructions in Makefile. Before it **gfortran** compiler must be installed on your computer! In case of successful compilation, an executable file **flex.exe** must appear on the directory above. To run the calculation program, in the same directory, along with **flex.exe**, there must be **INPUT.txt**-file, containing the parameters of the simulated system and instructions for executing the code. Its example is uploaded to the same directory as the **Source** folder.

## Input files
### INPUT.txt (required file)
This file contains following instructions:

| instructions | assignments |
|------------|-------------------------------------------------|
|**N1**      |  polymerization degree of chain A               |
|**sigma1**  |  grafting density of chain A                    |
|**p1**      |  Kuhn segment length of chain A                 |
|**chi1**    |  Flory's parameter for chain A / solvent        |
|**zp1**     |  z-coordinate of the grafting point (A-chain)   |
|**N2**      |  polymerization degree of chain B               |
|**sigma2**  |  grafting density of chain B                    |
|**p2**      |  Kuhn segment length of chain B                 |
|**chi2**    |  Flory's parameter for chain B / solvent        |
|**zp2**     |  z-coordinate of the grafting point (B-chain)   |
|**n_layer** |  number of layers                               |
|**left_w**  |  type of left wall (0 - solid, 1 - mirror)      |
|**right_w** |  type of right wall (0 - solid, 1 - mirror)     |
|**chi12**   |  Flory's parameter for polymer-polymer          |
|**eta**     |  step size of gradient descent for $\alpha(z)$  |
|**ksi**     |  step size of gradient descent for $u_{int}(z)$ |
|**nfree**   |  number of "free steps" at descent              |
|**swpro**   |  if swpro=0: switch off print of profiles       |

The order of listing these parameters in the **INPUT.txt**-file must be strictly observed without gaps between lines

### initial_guess.in (optional file)
**initial_guess.in** is both an input and an output file.

This file is automatically created/overwritten on every successful execution of the program. It contains information on the optimal convergence step, Lagrange field and chemical potential fields for polymer chains and solvent.

If this file is present in the directory along with **flex.exe**, the new calculation will be performed with the initial conditions specified in this file. This can significantly reduce the calculation time for large and complex systems. However, it can lead to incorrect results if the system is metastable and has several local minima of the free energy close in depth. Therefore, it is recommended to run the calculation without this file by default.

## Output files

### INFOR.info
Service information about the execution of the calculation algorithm is recorded to this file. Informs about the successful completion of the program, or gives error messages.

### data.out
Here the first moment values of both chains (H1 and H2) and free energy of system (F) are printing.

### profile.out
Here the volume fraction (phi1 and phi2) and free ends (end1 and end2) distributions of both chains are printing.


