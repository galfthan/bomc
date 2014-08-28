Bomc
====

BOMC is a Monte Carlo based generator of defect-free amorphous samples
of Si and SiO2, written by Sebastian von Alfthan (galfthan at iki dot fi).


If you publish research where Bomc has been utilized then you are
kindly asked to cite: "Realistic models of amorphous silica: A
comparative study of different potentials", S. von Alfthan,
A. Kuronen, and K. Kaski, Phys. Rev. B 68, 073203
(2003). ([pdf](http://lib.tkk.fi/Diss/2006/isbn9512285401/article2.pdf))



The code was developed as a personal research code in 2002 - 2003, and
thus there is little documentation. It is best to look at the code,
and the two example cases, to understand how it operates. The
[article](http://lib.tkk.fi/Diss/2006/isbn9512285401/article2.pdf)
mentioned above, and my
[thesis](http://lib.tkk.fi/Diss/2006/isbn9512285401/isbn9512285401.pdf)
also describes in some details the method. A more detailed discussion
of the optimized algorithms implemented in the code is available from
epaps: [Optimization of the
WWW-method](http://netserver.aip.org/epaps/phys_rev_b/E-PRBMDO-68-090327/epaps.pdf)

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

