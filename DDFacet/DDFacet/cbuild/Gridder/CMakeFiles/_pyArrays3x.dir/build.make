# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.18

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Disable VCS-based implicit rules.
% : %,v


# Disable VCS-based implicit rules.
% : RCS/%


# Disable VCS-based implicit rules.
% : RCS/%,v


# Disable VCS-based implicit rules.
% : SCCS/s.%


# Disable VCS-based implicit rules.
% : s.%


.SUFFIXES: .hpux_make_needs_suffix_list


# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /usr/local/src/DDFacet/DDFacet

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /usr/local/src/DDFacet/DDFacet/cbuild

# Include any dependencies generated for this target.
include Gridder/CMakeFiles/_pyArrays3x.dir/depend.make

# Include the progress variables for this target.
include Gridder/CMakeFiles/_pyArrays3x.dir/progress.make

# Include the compile flags for this target's objects.
include Gridder/CMakeFiles/_pyArrays3x.dir/flags.make

Gridder/CMakeFiles/_pyArrays3x.dir/Arrays.cc.o: Gridder/CMakeFiles/_pyArrays3x.dir/flags.make
Gridder/CMakeFiles/_pyArrays3x.dir/Arrays.cc.o: ../Gridder/Arrays.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/usr/local/src/DDFacet/DDFacet/cbuild/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object Gridder/CMakeFiles/_pyArrays3x.dir/Arrays.cc.o"
	cd /usr/local/src/DDFacet/DDFacet/cbuild/Gridder && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/_pyArrays3x.dir/Arrays.cc.o -c /usr/local/src/DDFacet/DDFacet/Gridder/Arrays.cc

Gridder/CMakeFiles/_pyArrays3x.dir/Arrays.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/_pyArrays3x.dir/Arrays.cc.i"
	cd /usr/local/src/DDFacet/DDFacet/cbuild/Gridder && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /usr/local/src/DDFacet/DDFacet/Gridder/Arrays.cc > CMakeFiles/_pyArrays3x.dir/Arrays.cc.i

Gridder/CMakeFiles/_pyArrays3x.dir/Arrays.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/_pyArrays3x.dir/Arrays.cc.s"
	cd /usr/local/src/DDFacet/DDFacet/cbuild/Gridder && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /usr/local/src/DDFacet/DDFacet/Gridder/Arrays.cc -o CMakeFiles/_pyArrays3x.dir/Arrays.cc.s

# Object files for target _pyArrays3x
_pyArrays3x_OBJECTS = \
"CMakeFiles/_pyArrays3x.dir/Arrays.cc.o"

# External object files for target _pyArrays3x
_pyArrays3x_EXTERNAL_OBJECTS =

Gridder/_pyArrays3x.so: Gridder/CMakeFiles/_pyArrays3x.dir/Arrays.cc.o
Gridder/_pyArrays3x.so: Gridder/CMakeFiles/_pyArrays3x.dir/build.make
Gridder/_pyArrays3x.so: /usr/lib/x86_64-linux-gnu/librt.so
Gridder/_pyArrays3x.so: /usr/lib/x86_64-linux-gnu/libpython3.9.so
Gridder/_pyArrays3x.so: Gridder/CMakeFiles/_pyArrays3x.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/usr/local/src/DDFacet/DDFacet/cbuild/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX shared library _pyArrays3x.so"
	cd /usr/local/src/DDFacet/DDFacet/cbuild/Gridder && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/_pyArrays3x.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
Gridder/CMakeFiles/_pyArrays3x.dir/build: Gridder/_pyArrays3x.so

.PHONY : Gridder/CMakeFiles/_pyArrays3x.dir/build

Gridder/CMakeFiles/_pyArrays3x.dir/clean:
	cd /usr/local/src/DDFacet/DDFacet/cbuild/Gridder && $(CMAKE_COMMAND) -P CMakeFiles/_pyArrays3x.dir/cmake_clean.cmake
.PHONY : Gridder/CMakeFiles/_pyArrays3x.dir/clean

Gridder/CMakeFiles/_pyArrays3x.dir/depend:
	cd /usr/local/src/DDFacet/DDFacet/cbuild && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /usr/local/src/DDFacet/DDFacet /usr/local/src/DDFacet/DDFacet/Gridder /usr/local/src/DDFacet/DDFacet/cbuild /usr/local/src/DDFacet/DDFacet/cbuild/Gridder /usr/local/src/DDFacet/DDFacet/cbuild/Gridder/CMakeFiles/_pyArrays3x.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : Gridder/CMakeFiles/_pyArrays3x.dir/depend

