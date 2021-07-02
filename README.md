# Global-TEP
Global-TEP is a global solver for TEP (transmission expansion planning) problem considering the AC model of the transmission network. Global-TEP codes are written in Julia 1.6 using JuMP v0.21.6 as the mathematical optimization modeling language. To solve the LB and UB models, ILOG CPLEX v20.1 and Ipopt v3.12 are respectively used.
To understand how Global-TEP works, you can refer to the manuscript: M. Mehrtash and Y. Cao, "A New Global Solver for Transmission Expansion Planning with AC Network Model," in IEEE Transactions on Power Systems, doi: 10.1109/TPWRS.2021.3086085.

## Developers
Global-TEP has been developed by:
    
    Mahdi Mehrtash, University of British Columbia
    Yankai Cao, University of British Columbia
    
## Citing Global-TEP
If you find Global-TEP useful for your work, you might cite the manuscript:

    @article{Mehrtash2021,
    author="Mehrtash, Mahdi
    and Cao, Yankai",
    title="A New Global Solver for Transmission Expansion Planning with AC Network Model",
    journal="IEEE Transactions on Power Systems",
    year="2021",
    month="June",
    day="2",
    issn="0885-8950",
    doi="10.1109/TPWRS.2021.3086085",
    url="https://ieeexplore.ieee.org/document/9445630"
    }

## Features
    compatible with Julia 1.6
    compatible with the lastest version of JuMP
    Support MINLP
    Support parallel computing
    
