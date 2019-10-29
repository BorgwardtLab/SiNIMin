CMAKE_MINIMUM_REQUIRED( VERSION 3.0 )
PROJECT( SiNIMin CXX )

SET( CMAKE_CXX_STANDARD 11 )

FIND_PACKAGE( Eigen3 REQUIRED )
FIND_PACKAGE( OpenMP COMPONENTS CXX REQUIRED )

ADD_EXECUTABLE( sinimin
  main/sinimin.cpp
)

ADD_EXECUTABLE( sinimin_fwer
  main/sinimin_fwer_estimation.cpp
)

TARGET_COMPILE_OPTIONS( sinimin      PUBLIC "-O5" )
TARGET_COMPILE_OPTIONS( sinimin_fwer PUBLIC "-O5" )

TARGET_LINK_LIBRARIES( sinimin      Eigen3::Eigen OpenMP::OpenMP_CXX )
TARGET_LINK_LIBRARIES( sinimin_fwer Eigen3::Eigen OpenMP::OpenMP_CXX )

INSTALL( TARGETS sinimin sinimin_fwer
  RUNTIME
  DESTINATION bin
)
