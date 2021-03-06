# MBT-Coffee
MBT-Coffee, is based on the MEL library reduction and incorporates a new objective function based on the use of both the consistency library and a substitution matrix in order to improve the final alignment accuracy and also providing a great reduction in the execution time.

Prerequisites
--------------
T-Coffee compilation requires the following tools installed on your system ``make`` and ``gcc-c++``. 


Compile 
--------
Download the git repository on your computer.
    
Make sure you have installed the required dependencies listed above. 
When done, move in the project root folder and enter the following commands:     
    
    $ cd src
    $ make
    

The binary will be automatically generated in the folder.


Example
--------

It is included a script named ``run_example`` which executes t_coffee with the required parameters.
The example uses the ``GEL_100`` dataset from ``HomFam``.

In order to use the MEL reduction, download ``TCoffee-MEL`` or ``PPCAS-MEL`` to generate a reduced library (-out_lib parameter with TCoffee-MEL) and use it as input in ``MBT-Coffee`` (-lib parameter).

