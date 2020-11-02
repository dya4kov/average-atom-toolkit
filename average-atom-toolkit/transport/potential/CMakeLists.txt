set(POTENTIAL_SOURCES
  ${CMAKE_CURRENT_SOURCE_DIR}/base.cxx
  ${CMAKE_CURRENT_SOURCE_DIR}/debye.cxx
  ${CMAKE_CURRENT_SOURCE_DIR}/hard-sphere.cxx
  ${CMAKE_CURRENT_SOURCE_DIR}/polarization.cxx
  ${CMAKE_CURRENT_SOURCE_DIR}/square-well.cxx
  ${CMAKE_CURRENT_SOURCE_DIR}/tabular.cxx
    PARENT_SCOPE
)
# install(FILES   base.h
# 	            debye.h
# 	            hard-sphere.h
# 	            polarization.h
# 	            square-well.h
# 	            tabular.h
#     DESTINATION include/${PROJECT_NAME}/${LIB_DIR_NAME}/potential
# )