cmake_minimum_required(VERSION 2.6)
find_package(Boost 1.36.0)
find_library(CONFIG++_LIBRARY NAMES config++ PATH /usr/lib /usr/local/lib)

if(Boost_FOUND)
  include_directories(${Boost_INCLUDE_DIRS})
  include_directories(''/usr/local/include'')
  add_executable(Polymer temp_pf.cpp)
  target_link_libraries(Polymer ${CONFIG++_LIBRARY})
endif()