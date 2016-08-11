FILE(GLOB headers "${CMAKE_CURRENT_SOURCE_DIR}/*.hpp")
FILE(GLOB sources "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp")

ADD_EXECUTABLE(liegroups_test_expm expm.cpp)
ADD_EXECUTABLE(liegroups_test_extrinsics extrinsics.cpp)
ADD_EXECUTABLE(liegroups_test_main main.cpp)
ADD_EXECUTABLE(liegroups_test_moments moments.cpp)
ADD_EXECUTABLE(liegroups_test_se3_gen se3_gen.cpp)
ADD_EXECUTABLE(liegroups_test_se3_moments se3_moments.cpp)
ADD_EXECUTABLE(liegroups_test_spline_test spline_test.cpp)

TARGET_LINK_LIBRARIES(liegroups_test_expm liegroups)
TARGET_LINK_LIBRARIES(liegroups_test_extrinsics liegroups)
TARGET_LINK_LIBRARIES(liegroups_test_main liegroups)
TARGET_LINK_LIBRARIES(liegroups_test_moments liegroups)
TARGET_LINK_LIBRARIES(liegroups_test_se3_gen liegroups)
TARGET_LINK_LIBRARIES(liegroups_test_se3_moments liegroups)
TARGET_LINK_LIBRARIES(liegroups_test_spline_test liegroups)

INSTALL(TARGETS
        liegroups_test_expm
        liegroups_test_extrinsics
        liegroups_test_main
        liegroups_test_moments
        liegroups_test_se3_gen
        liegroups_test_se3_moments
        liegroups_test_spline_test
  RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/bin
  LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
  ARCHIVE DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
)

INSTALL(FILES ${headers} DESTINATION ${CMAKE_INSTALL_PREFIX}/include/liegroups/test)
