# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.25

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
CMAKE_COMMAND = /home/gaxirs/anaconda3/lib/python3.9/site-packages/cmake/data/bin/cmake

# The command to remove a file.
RM = /home/gaxirs/anaconda3/lib/python3.9/site-packages/cmake/data/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = "/home/gaxirs/Documents/Cours/Q6/LEPL1110/Projets/Project with graphs"

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = "/home/gaxirs/Documents/Cours/Q6/LEPL1110/Projets/Project with graphs/build"

# Include any dependencies generated for this target.
include CMakeFiles/myFem.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/myFem.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/myFem.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/myFem.dir/flags.make

CMakeFiles/myFem.dir/src/fem.c.o: CMakeFiles/myFem.dir/flags.make
CMakeFiles/myFem.dir/src/fem.c.o: /home/gaxirs/Documents/Cours/Q6/LEPL1110/Projets/Project\ with\ graphs/src/fem.c
CMakeFiles/myFem.dir/src/fem.c.o: CMakeFiles/myFem.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir="/home/gaxirs/Documents/Cours/Q6/LEPL1110/Projets/Project with graphs/build/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_1) "Building C object CMakeFiles/myFem.dir/src/fem.c.o"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -MD -MT CMakeFiles/myFem.dir/src/fem.c.o -MF CMakeFiles/myFem.dir/src/fem.c.o.d -o CMakeFiles/myFem.dir/src/fem.c.o -c "/home/gaxirs/Documents/Cours/Q6/LEPL1110/Projets/Project with graphs/src/fem.c"

CMakeFiles/myFem.dir/src/fem.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/myFem.dir/src/fem.c.i"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E "/home/gaxirs/Documents/Cours/Q6/LEPL1110/Projets/Project with graphs/src/fem.c" > CMakeFiles/myFem.dir/src/fem.c.i

CMakeFiles/myFem.dir/src/fem.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/myFem.dir/src/fem.c.s"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S "/home/gaxirs/Documents/Cours/Q6/LEPL1110/Projets/Project with graphs/src/fem.c" -o CMakeFiles/myFem.dir/src/fem.c.s

CMakeFiles/myFem.dir/src/glfem.c.o: CMakeFiles/myFem.dir/flags.make
CMakeFiles/myFem.dir/src/glfem.c.o: /home/gaxirs/Documents/Cours/Q6/LEPL1110/Projets/Project\ with\ graphs/src/glfem.c
CMakeFiles/myFem.dir/src/glfem.c.o: CMakeFiles/myFem.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir="/home/gaxirs/Documents/Cours/Q6/LEPL1110/Projets/Project with graphs/build/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_2) "Building C object CMakeFiles/myFem.dir/src/glfem.c.o"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -MD -MT CMakeFiles/myFem.dir/src/glfem.c.o -MF CMakeFiles/myFem.dir/src/glfem.c.o.d -o CMakeFiles/myFem.dir/src/glfem.c.o -c "/home/gaxirs/Documents/Cours/Q6/LEPL1110/Projets/Project with graphs/src/glfem.c"

CMakeFiles/myFem.dir/src/glfem.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/myFem.dir/src/glfem.c.i"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E "/home/gaxirs/Documents/Cours/Q6/LEPL1110/Projets/Project with graphs/src/glfem.c" > CMakeFiles/myFem.dir/src/glfem.c.i

CMakeFiles/myFem.dir/src/glfem.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/myFem.dir/src/glfem.c.s"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S "/home/gaxirs/Documents/Cours/Q6/LEPL1110/Projets/Project with graphs/src/glfem.c" -o CMakeFiles/myFem.dir/src/glfem.c.s

CMakeFiles/myFem.dir/src/solve.c.o: CMakeFiles/myFem.dir/flags.make
CMakeFiles/myFem.dir/src/solve.c.o: /home/gaxirs/Documents/Cours/Q6/LEPL1110/Projets/Project\ with\ graphs/src/solve.c
CMakeFiles/myFem.dir/src/solve.c.o: CMakeFiles/myFem.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir="/home/gaxirs/Documents/Cours/Q6/LEPL1110/Projets/Project with graphs/build/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_3) "Building C object CMakeFiles/myFem.dir/src/solve.c.o"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -MD -MT CMakeFiles/myFem.dir/src/solve.c.o -MF CMakeFiles/myFem.dir/src/solve.c.o.d -o CMakeFiles/myFem.dir/src/solve.c.o -c "/home/gaxirs/Documents/Cours/Q6/LEPL1110/Projets/Project with graphs/src/solve.c"

CMakeFiles/myFem.dir/src/solve.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/myFem.dir/src/solve.c.i"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E "/home/gaxirs/Documents/Cours/Q6/LEPL1110/Projets/Project with graphs/src/solve.c" > CMakeFiles/myFem.dir/src/solve.c.i

CMakeFiles/myFem.dir/src/solve.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/myFem.dir/src/solve.c.s"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S "/home/gaxirs/Documents/Cours/Q6/LEPL1110/Projets/Project with graphs/src/solve.c" -o CMakeFiles/myFem.dir/src/solve.c.s

CMakeFiles/myFem.dir/src/main.c.o: CMakeFiles/myFem.dir/flags.make
CMakeFiles/myFem.dir/src/main.c.o: /home/gaxirs/Documents/Cours/Q6/LEPL1110/Projets/Project\ with\ graphs/src/main.c
CMakeFiles/myFem.dir/src/main.c.o: CMakeFiles/myFem.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir="/home/gaxirs/Documents/Cours/Q6/LEPL1110/Projets/Project with graphs/build/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_4) "Building C object CMakeFiles/myFem.dir/src/main.c.o"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -MD -MT CMakeFiles/myFem.dir/src/main.c.o -MF CMakeFiles/myFem.dir/src/main.c.o.d -o CMakeFiles/myFem.dir/src/main.c.o -c "/home/gaxirs/Documents/Cours/Q6/LEPL1110/Projets/Project with graphs/src/main.c"

CMakeFiles/myFem.dir/src/main.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/myFem.dir/src/main.c.i"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E "/home/gaxirs/Documents/Cours/Q6/LEPL1110/Projets/Project with graphs/src/main.c" > CMakeFiles/myFem.dir/src/main.c.i

CMakeFiles/myFem.dir/src/main.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/myFem.dir/src/main.c.s"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S "/home/gaxirs/Documents/Cours/Q6/LEPL1110/Projets/Project with graphs/src/main.c" -o CMakeFiles/myFem.dir/src/main.c.s

# Object files for target myFem
myFem_OBJECTS = \
"CMakeFiles/myFem.dir/src/fem.c.o" \
"CMakeFiles/myFem.dir/src/glfem.c.o" \
"CMakeFiles/myFem.dir/src/solve.c.o" \
"CMakeFiles/myFem.dir/src/main.c.o"

# External object files for target myFem
myFem_EXTERNAL_OBJECTS =

myFem: CMakeFiles/myFem.dir/src/fem.c.o
myFem: CMakeFiles/myFem.dir/src/glfem.c.o
myFem: CMakeFiles/myFem.dir/src/solve.c.o
myFem: CMakeFiles/myFem.dir/src/main.c.o
myFem: CMakeFiles/myFem.dir/build.make
myFem: glfw/src/libglfw3.a
myFem: /usr/lib/x86_64-linux-gnu/libGL.so
myFem: /home/gaxirs/Documents/Cours/Q6/LEPL1110/Projets/Project\ with\ graphs/gmsh/gmsh-4.11.1-Linux64-sdk/lib/libgmsh.so
myFem: /usr/lib/x86_64-linux-gnu/librt.a
myFem: /usr/lib/x86_64-linux-gnu/libm.so
myFem: /usr/lib/x86_64-linux-gnu/libX11.so
myFem: CMakeFiles/myFem.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir="/home/gaxirs/Documents/Cours/Q6/LEPL1110/Projets/Project with graphs/build/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_5) "Linking C executable myFem"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/myFem.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/myFem.dir/build: myFem
.PHONY : CMakeFiles/myFem.dir/build

CMakeFiles/myFem.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/myFem.dir/cmake_clean.cmake
.PHONY : CMakeFiles/myFem.dir/clean

CMakeFiles/myFem.dir/depend:
	cd "/home/gaxirs/Documents/Cours/Q6/LEPL1110/Projets/Project with graphs/build" && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" "/home/gaxirs/Documents/Cours/Q6/LEPL1110/Projets/Project with graphs" "/home/gaxirs/Documents/Cours/Q6/LEPL1110/Projets/Project with graphs" "/home/gaxirs/Documents/Cours/Q6/LEPL1110/Projets/Project with graphs/build" "/home/gaxirs/Documents/Cours/Q6/LEPL1110/Projets/Project with graphs/build" "/home/gaxirs/Documents/Cours/Q6/LEPL1110/Projets/Project with graphs/build/CMakeFiles/myFem.dir/DependInfo.cmake" --color=$(COLOR)
.PHONY : CMakeFiles/myFem.dir/depend

