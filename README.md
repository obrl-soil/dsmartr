# dsmartr
An R implementation of the DSMART algorithm

This package builds on [some earlier work](https://github.com/obrl-soil/disaggregation) adapting 
[the original R-based version of dsmart](https://bitbucket.org/brendo1001/dsmart) 
for my own needs. This implementation incorporates extensive changes to data handling and sampling and relies on the sf package over sp as
far as is it can. 

Install with  

    devtools::install_github("obrl-soil/dsmartr")

NB: Since sf is under active development, install will force an upgrade of sf to >= 0.5.
