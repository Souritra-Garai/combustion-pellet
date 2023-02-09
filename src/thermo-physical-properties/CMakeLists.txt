add_library(CoreShellParticle Core-Shell-Particle.cpp)
add_library(ThermalConductivityPellet Thermal-Conductivity-Pellet.cpp)

add_library(PackedPellet Packed-Pellet.cpp)
target_link_libraries(PackedPellet PUBLIC ThermalConductivityPellet CoreShellParticle ReadData)