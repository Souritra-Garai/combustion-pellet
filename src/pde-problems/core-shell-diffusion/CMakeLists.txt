add_library(CoreShellDiffusion Utilities.cpp Equations.cpp)
target_link_libraries(CoreShellDiffusion PUBLIC CoreShellParticle LUSolver ReadData)