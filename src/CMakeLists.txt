# CMakeList.txt : CMake project for CMakeProject1, include source and define
# project specific logic here.
#
cmake_minimum_required (VERSION 3.10)

## QT5
find_package(Qt5 COMPONENTS Core Gui Widgets REQUIRED)
message("Qt path: " ${Qt5_DIR})
set_property(GLOBAL PROPERTY USE_FOLDERS ON)
set_property(GLOBAL PROPERTY AUTOGEN_TARGETS_FOLDER "${CMAKE_SOURCE_DIR}")
set (CMAKE_AUTOMOC ON)
set (CMAKE_AUTORCC ON)
set (CMAKE_AUTOUIC ON)
if(CMAKE_VERSION VERSION_LESS "3.7.0")
	set(CMAKE_INCLUDE_CURRENT_DIR ON)
endif()
set (QT5MODULE Qt5::Core Qt5::Gui Qt5::Widgets)
#qt5_wrap_ui(WRAP_FILES ${file_ui})
if (POLICY CMP0100)
	cmake_policy(SET CMP0100 NEW)
endif()
# include ui source files
FILE(GLOB_RECURSE file_ui ./ui/*.ui)
# message("ui files: " ${file_ui})

## VTK
find_package(VTK REQUIRED)
message ("VTK found. VERSION: " ${VTK_VERSION})
if (VTK_VERSION VERSION_LESS "8.90.0")
	include (${VTK_USE_FILE})
endif()

## Eigen3
if (NOT Eigen3_FOUND)
	if (UNIX)
		find_package (Eigen3 REQUIRED)
		include_directories(${EIGEN3_INCLUDE_DIR})
		message(${EIGEN3_INCLUDE_DIR})
	endif()
	if(WIN32)
		include_directories("$ENV{EIGEN_DIR}")
		message("Eigen3 directory : " "$ENV{EIGEN_DIR}")
		set(Eigen3_FOUND true)
	endif()
endif()


## source files
if(${CMAKE_VERSION} VERSION_GREATER 3.12)
	file(GLOB_RECURSE src_files *.cc *.cxx *.cpp *.h *.hh *.hpp CONFIG_DEPENDS)
else()
	file(GLOB_RECURSE src_files *.cc *.cxx *.cpp *.h *.hh *.hpp)
endif()
# file(GLOB_RECURSE OVM_SRC ./OpenVolumeMesh "*.hh" "*.cc")
set(source_list ${src_files} ${file_ui})

add_executable(${program_name} ${source_list})

target_include_directories(${program_name} 
PUBLIC 
	${CMAKE_CURRENT_SOURCE_DIR}
	${EIGEN3_INCLUDE_DIR}
)
target_compile_features(${program_name} PRIVATE cxx_std_14)
set_target_properties(${program_name} PROPERTIES CXX_EXTENSIONS OFF)
## OpenVolumeMesh
if (MSVC)
## gcc, clang has no bigobj problem
	add_compile_options(/bigobj)
endif()
## define INCLUDE_TEMPLATES to enable OVM templates functions, like FileManager
target_compile_definitions(${program_name} PRIVATE -DINCLUDE_TEMPLATES)

if(MSVC)
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE ${PROJECT_SOURCE_DIR}/../out/bin) 
endif()

## add link libraries
target_link_libraries(${program_name} PRIVATE
    ${QT5MODULE}
    ${VTK_LIBRARIES}
)

# vtk autoinit
if (VTK_VERSION VERSION_GREATER_EQUAL "8.90.0")
	# vtk_module_autoinit is needed	
	vtk_module_autoinit(
		TARGETS ${program_name}
		MODULES ${VTK_LIBRARIES}
	)
endif()

# group source files in Visual Studio in folder way
if(MSVC)
source_group(TREE ${CMAKE_CURRENT_SOURCE_DIR} FILES ${source_list})
# foreach(srcf ${source_list})
	# file(RELATIVE_PATH source "${PROJECT_SOURCE_DIR}" "${srcf}")
	# get_filename_component(source_path "${source}" DIRECTORY)
	# string(REPLACE "/" "\\" source_path_msvc "${source_path}")
	# source_group(${source_path_msvc} FILES ${source})
# endforeach()
endif()

# Add tests and install targets if needed.  