Bomc
====

BOMC is a Monte Carlo based generator of defect-free amorphous samples
of Si and SiO2. 

If you publish research where Bomc has been utilized then you are
kindly asked to cite: S. von Alfthan, A. Kuronen, and K. Kaski,
Phys. Rev. B 68, 073203 (2003).


The code was developed as a research code in 2002 - 2003 for personal
use, and thus there is little documentation. It is best to look at the
code, and the two example cases, to understand how it operates.


## Compiling:

On a standard unix machine the code should compile by running make in
the src folder. One can edit the Makefile if it does not
compile. 

There is also a naive and not-very-well-tested parallel version that
one can enable in the makefile; scalability is probably limited to a
node or so.


## Usage:

When you run the program you should have three .par files in your
folder where all settings are defined:

* *main.par* General settings, the testcases have some comments in
them which should explain most options 
* *cg.par* Options for the optimization algorithms. There is a number of heuristics, some of
which are explained in the article. One probably do not need to touch
this too much. 
* *wwwPot.par* Potential model parameters

When the code is run, no parameters are given on the command line:

`./bomc`

