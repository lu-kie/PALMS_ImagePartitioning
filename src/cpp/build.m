% Build mex
 mex CXXFLAGS='$CXXFLAGS -fopenmp' LDFLAGS='-larmadillo -fopenmp'...
 	 LinewiseSolver_mexWrapper.cpp ArmadilloConverter.cpp...
     Compute1rErrors.cpp Extract1Dstripes.cpp FindBest1DPartition.cpp GenerateSystemMatrices.cpp...
     GetIndexes.cpp Interval.cpp LinewisePartitioning.cpp ReconstructionFromPartition.cpp Stripe.cpp
