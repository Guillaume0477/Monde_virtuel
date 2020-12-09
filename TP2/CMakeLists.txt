cmake_minimum_required(VERSION 3.7.0)

project(pgm LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTOUIC ON)


find_package(Qt5 COMPONENTS Widgets OpenGL REQUIRED)
find_package(OpenGL REQUIRED)

find_package(GLEW)
if(NOT ${GLEW_FOUND})

  message(WARNING "Glew not found, so manual compilation of glew added to the project")
  
  option(glew-cmake_BUILD_SHARED "Build the shared glew library" OFF)
  option(ONLY_LIBS "Do not build executables" ON)

  add_subdirectory(extern/glew)
  include_directories(${PROJECT_PATH}/extern/glew/include)
  set(USED_GLEW_LIBRARY libglew_static)
else()
  include_directories(${GLEW_INCLUDE_DIRS})
  set(USED_GLEW_LIBRARY ${GLEW_LIBRARIES})
endif()

file(
GLOB_RECURSE
source_files
"main.cpp"
)

add_executable( 
  pgm
  ${source_files}
)

target_link_libraries(pgm ${USED_GLEW_LIBRARY} ${OPENGL_LIBRARIES} Qt5::Widgets Qt5::OpenGL	)
