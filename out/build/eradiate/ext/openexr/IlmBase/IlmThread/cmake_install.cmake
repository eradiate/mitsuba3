# Install script for directory: /home/matt/rayference/eradiate/ext/mitsuba/ext/openexr/IlmBase/IlmThread

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/usr/local")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "Release")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "1")
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

# Set default install directory permissions.
if(NOT DEFINED CMAKE_OBJDUMP)
  set(CMAKE_OBJDUMP "/usr/bin/llvm-objdump-11")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libIlmThread-mitsuba.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libIlmThread-mitsuba.so")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libIlmThread-mitsuba.so"
         RPATH "$ORIGIN:$ORIGIN/../drjit")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE SHARED_LIBRARY FILES "/home/matt/rayference/eradiate/ext/mitsuba/out/build/eradiate/libIlmThread-mitsuba.so")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libIlmThread-mitsuba.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libIlmThread-mitsuba.so")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libIlmThread-mitsuba.so"
         OLD_RPATH "/home/matt/rayference/eradiate/ext/mitsuba/out/build/eradiate:"
         NEW_RPATH "$ORIGIN:$ORIGIN/../drjit")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/llvm-strip-11" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libIlmThread-mitsuba.so")
    endif()
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/OpenEXR" TYPE FILE FILES
    "/home/matt/rayference/eradiate/ext/mitsuba/ext/openexr/IlmBase/IlmThread/IlmThreadPool.h"
    "/home/matt/rayference/eradiate/ext/mitsuba/ext/openexr/IlmBase/IlmThread/IlmThread.h"
    "/home/matt/rayference/eradiate/ext/mitsuba/ext/openexr/IlmBase/IlmThread/IlmThreadSemaphore.h"
    "/home/matt/rayference/eradiate/ext/mitsuba/ext/openexr/IlmBase/IlmThread/IlmThreadMutex.h"
    "/home/matt/rayference/eradiate/ext/mitsuba/ext/openexr/IlmBase/IlmThread/IlmThreadNamespace.h"
    "/home/matt/rayference/eradiate/ext/mitsuba/ext/openexr/IlmBase/IlmThread/IlmThreadExport.h"
    "/home/matt/rayference/eradiate/ext/mitsuba/ext/openexr/IlmBase/IlmThread/IlmThreadForward.h"
    )
endif()

