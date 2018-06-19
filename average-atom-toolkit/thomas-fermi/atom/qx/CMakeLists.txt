add_library(ThomasFermiAtomQX SHARED
    potential.cxx
)
target_link_libraries(
  ThomasFermiAtomQX
  ThomasFermiChemicalPotential
  ThomasFermiChemicalPotentialQX
  ${NUM_TOOLS_LIBS}
)
install(TARGETS ThomasFermiAtomQX DESTINATION lib)
install(FILES potential.h
        DESTINATION include/average-atom-toolkit/thomas-fermi/atom/qx
)