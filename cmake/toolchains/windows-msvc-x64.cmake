set(CMAKE_SYSTEM_NAME "Windows" CACHE STRING "System name is Windows")
set(CMAKE_SYSREM_PROCESSOR "x86_64" CACHE STRING "System processor is x86_64")

set(GIT_EXECUTABLE "C:/Program Files/Microsoft Visual Studio/2022/Community/Common7/IDE/CommonExtensions/Microsoft/TeamFoundation/Team Explorer/Git/cmd/git.exe" CACHE FILEPATH "Path to git executable")

set(CMAKE_C_COMPILER_LAUNCHER "ccache.exe" CACHE FILEPATH "Enable ccache")
set(CMAKE_CXX_COMPILER_LAUNCHER "ccache.exe" CACHE FILEPATH "Enable ccache")

set(CMAKE_C_COMPILER "cl.exe" CACHE FILEPATH "MSVC C compiler")
set(CMAKE_CXX_COMPILER "cl.exe" CACHE FILEPATH "MSVC C++ compiler")

set(CMAKE_LINKER_TYPE "MSVC" CACHE STRING "MSVC linker")

set(CMAKE_C_FLAGS_INIT "/arch:AVX2 /EHsc /utf-8 /permissive-" CACHE STRING "C compiler flags")
set(CMAKE_CXX_FLAGS_INIT "${CMAKE_C_FLAGS_INIT}" CACHE STRING "C++ compiler flags")

set(CMAKE_C_FLAGS_DEBUG "/Od /Zi -D_DEBUG" CACHE STRING "C compiler flags for debug builds")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG}" CACHE STRING "C++ compiler flags for debug builds")
set(CMAKE_INTERPROCEDURAL_OPTIMIZATION_DEBUG False CACHE BOOL "Interprocedural optimization for debug builds")

set(CMAKE_C_FLAGS_RELEASE "/O2" CACHE STRING "C compiler flags for release builds")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE}" CACHE STRING "C++ compiler flags for release builds")
set(CMAKE_INTERPROCEDURAL_OPTIMIZATION_RELEASE True CACHE BOOL "Interprocedural optimization for release builds")

set(CMAKE_C_FLAGS_RELWITHDEBINFO "/O2 /Zi" CACHE STRING "C compiler flags for release with debug info builds")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_C_FLAGS_RELWITHDEBINFO}" CACHE STRING "C++ compiler flags for release with debug info builds")
set(CMAKE_INTERPROCEDURAL_OPTIMIZATION_RELWITHDEBINFO True CACHE BOOL "Interprocedural optimization for release with debug info builds")

set(CMAKE_C_FLAGS_MINSIZEREL "/Os" CACHE STRING "C compiler flags for minimum size release builds")
set(CMAKE_CXX_FLAGS_MINSIZEREL "${CMAKE_C_FLAGS_MINSIZEREL}" CACHE STRING "C++ compiler flags for minimum size release builds")
set(CMAKE_INTERPROCEDURAL_OPTIMIZATION_MINSIZEREL True CACHE BOOL "Interprocedural optimization for minimum size release builds")

if (POLICY CMP0141)
  cmake_policy(SET CMP0141 NEW)
  set(CMAKE_MSVC_DEBUG_INFORMATION_FORMAT "$<IF:$<AND:$<C_COMPILER_ID:MSVC>,$<CXX_COMPILER_ID:MSVC>>,$<$<CONFIG:Debug,RelWithDebInfo>:EditAndContinue>,$<$<CONFIG:Debug,RelWithDebInfo>:ProgramDatabase>>" CACHE INTERNAL "Enable warm reload")
endif()
