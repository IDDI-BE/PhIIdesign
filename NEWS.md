## CHANGES IN PhIIdesign VERSION 0.2

- Added tinytest unit tests for 
    - fleming1stage
    - exact1stage
    - sargent1stage
    - sargent2stage
    - simon2stage
- Added website
- Added continuous integration using github workflows
- Makes sure R CMD check runs fine
- Added benchmarking script
- Make documentation consistent, put the commented examples which were in 0.1 code as examples in the functions + add example datasets of settings as datasets in the package
- Split up simon2stage in calculation and plotting
- Added speedier variant of 
    - fleming1stage, using Rcpp - speedup x17
        + addded fleming1stage_multiple which uses this speedup and make calling the function for several settings more easy
    - exact1stage, using Rcpp - speedup x40
        + added exact1stage_multiple which uses this speedup and make calling the function for several settings more easy
    - sargent1stage, using Rcpp - speedup 14% due to use of sargent1stage_N_r_s
        + sargent1stage uses this speedup and make calling the function for several settings more easy
    - simon2stage, using lookup - speedup x4
    - sargent2stage, using lookup similar as simon2stage - speedup x4
    - plots of probability of success e.g.
- Added allsinglearm

## CHANGES IN PhIIdesign VERSION 0.1

- Initial version of the package


