# Install script for directory: /home/matt/rayference/eradiate/ext/mitsuba/src/bsdfs

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
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/plugins/blendbsdf.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/plugins/blendbsdf.so")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/plugins/blendbsdf.so"
         RPATH "\$ORIGIN:\$ORIGIN/../drjit")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin/plugins" TYPE SHARED_LIBRARY FILES "/home/matt/rayference/eradiate/ext/mitsuba/out/build/eradiate/plugins/blendbsdf.so")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/plugins/blendbsdf.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/plugins/blendbsdf.so")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/plugins/blendbsdf.so"
         OLD_RPATH "/home/matt/rayference/eradiate/ext/mitsuba/out/build/eradiate:"
         NEW_RPATH "\$ORIGIN:\$ORIGIN/../drjit")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/llvm-strip-11" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/plugins/blendbsdf.so")
    endif()
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/plugins/bumpmap.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/plugins/bumpmap.so")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/plugins/bumpmap.so"
         RPATH "\$ORIGIN:\$ORIGIN/../drjit")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin/plugins" TYPE SHARED_LIBRARY FILES "/home/matt/rayference/eradiate/ext/mitsuba/out/build/eradiate/plugins/bumpmap.so")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/plugins/bumpmap.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/plugins/bumpmap.so")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/plugins/bumpmap.so"
         OLD_RPATH "/home/matt/rayference/eradiate/ext/mitsuba/out/build/eradiate:"
         NEW_RPATH "\$ORIGIN:\$ORIGIN/../drjit")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/llvm-strip-11" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/plugins/bumpmap.so")
    endif()
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/plugins/conductor.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/plugins/conductor.so")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/plugins/conductor.so"
         RPATH "\$ORIGIN:\$ORIGIN/../drjit")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin/plugins" TYPE SHARED_LIBRARY FILES "/home/matt/rayference/eradiate/ext/mitsuba/out/build/eradiate/plugins/conductor.so")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/plugins/conductor.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/plugins/conductor.so")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/plugins/conductor.so"
         OLD_RPATH "/home/matt/rayference/eradiate/ext/mitsuba/out/build/eradiate:"
         NEW_RPATH "\$ORIGIN:\$ORIGIN/../drjit")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/llvm-strip-11" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/plugins/conductor.so")
    endif()
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/plugins/dielectric.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/plugins/dielectric.so")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/plugins/dielectric.so"
         RPATH "\$ORIGIN:\$ORIGIN/../drjit")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin/plugins" TYPE SHARED_LIBRARY FILES "/home/matt/rayference/eradiate/ext/mitsuba/out/build/eradiate/plugins/dielectric.so")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/plugins/dielectric.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/plugins/dielectric.so")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/plugins/dielectric.so"
         OLD_RPATH "/home/matt/rayference/eradiate/ext/mitsuba/out/build/eradiate:"
         NEW_RPATH "\$ORIGIN:\$ORIGIN/../drjit")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/llvm-strip-11" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/plugins/dielectric.so")
    endif()
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/plugins/diffuse.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/plugins/diffuse.so")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/plugins/diffuse.so"
         RPATH "\$ORIGIN:\$ORIGIN/../drjit")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin/plugins" TYPE SHARED_LIBRARY FILES "/home/matt/rayference/eradiate/ext/mitsuba/out/build/eradiate/plugins/diffuse.so")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/plugins/diffuse.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/plugins/diffuse.so")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/plugins/diffuse.so"
         OLD_RPATH "/home/matt/rayference/eradiate/ext/mitsuba/out/build/eradiate:"
         NEW_RPATH "\$ORIGIN:\$ORIGIN/../drjit")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/llvm-strip-11" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/plugins/diffuse.so")
    endif()
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/plugins/hair.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/plugins/hair.so")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/plugins/hair.so"
         RPATH "\$ORIGIN:\$ORIGIN/../drjit")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin/plugins" TYPE SHARED_LIBRARY FILES "/home/matt/rayference/eradiate/ext/mitsuba/out/build/eradiate/plugins/hair.so")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/plugins/hair.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/plugins/hair.so")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/plugins/hair.so"
         OLD_RPATH "/home/matt/rayference/eradiate/ext/mitsuba/out/build/eradiate:"
         NEW_RPATH "\$ORIGIN:\$ORIGIN/../drjit")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/llvm-strip-11" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/plugins/hair.so")
    endif()
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/plugins/mask.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/plugins/mask.so")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/plugins/mask.so"
         RPATH "\$ORIGIN:\$ORIGIN/../drjit")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin/plugins" TYPE SHARED_LIBRARY FILES "/home/matt/rayference/eradiate/ext/mitsuba/out/build/eradiate/plugins/mask.so")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/plugins/mask.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/plugins/mask.so")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/plugins/mask.so"
         OLD_RPATH "/home/matt/rayference/eradiate/ext/mitsuba/out/build/eradiate:"
         NEW_RPATH "\$ORIGIN:\$ORIGIN/../drjit")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/llvm-strip-11" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/plugins/mask.so")
    endif()
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/plugins/measured.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/plugins/measured.so")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/plugins/measured.so"
         RPATH "\$ORIGIN:\$ORIGIN/../drjit")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin/plugins" TYPE SHARED_LIBRARY FILES "/home/matt/rayference/eradiate/ext/mitsuba/out/build/eradiate/plugins/measured.so")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/plugins/measured.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/plugins/measured.so")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/plugins/measured.so"
         OLD_RPATH "/home/matt/rayference/eradiate/ext/mitsuba/out/build/eradiate:"
         NEW_RPATH "\$ORIGIN:\$ORIGIN/../drjit")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/llvm-strip-11" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/plugins/measured.so")
    endif()
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/plugins/normalmap.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/plugins/normalmap.so")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/plugins/normalmap.so"
         RPATH "\$ORIGIN:\$ORIGIN/../drjit")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin/plugins" TYPE SHARED_LIBRARY FILES "/home/matt/rayference/eradiate/ext/mitsuba/out/build/eradiate/plugins/normalmap.so")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/plugins/normalmap.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/plugins/normalmap.so")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/plugins/normalmap.so"
         OLD_RPATH "/home/matt/rayference/eradiate/ext/mitsuba/out/build/eradiate:"
         NEW_RPATH "\$ORIGIN:\$ORIGIN/../drjit")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/llvm-strip-11" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/plugins/normalmap.so")
    endif()
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/plugins/null.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/plugins/null.so")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/plugins/null.so"
         RPATH "\$ORIGIN:\$ORIGIN/../drjit")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin/plugins" TYPE SHARED_LIBRARY FILES "/home/matt/rayference/eradiate/ext/mitsuba/out/build/eradiate/plugins/null.so")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/plugins/null.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/plugins/null.so")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/plugins/null.so"
         OLD_RPATH "/home/matt/rayference/eradiate/ext/mitsuba/out/build/eradiate:"
         NEW_RPATH "\$ORIGIN:\$ORIGIN/../drjit")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/llvm-strip-11" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/plugins/null.so")
    endif()
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/plugins/plastic.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/plugins/plastic.so")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/plugins/plastic.so"
         RPATH "\$ORIGIN:\$ORIGIN/../drjit")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin/plugins" TYPE SHARED_LIBRARY FILES "/home/matt/rayference/eradiate/ext/mitsuba/out/build/eradiate/plugins/plastic.so")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/plugins/plastic.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/plugins/plastic.so")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/plugins/plastic.so"
         OLD_RPATH "/home/matt/rayference/eradiate/ext/mitsuba/out/build/eradiate:"
         NEW_RPATH "\$ORIGIN:\$ORIGIN/../drjit")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/llvm-strip-11" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/plugins/plastic.so")
    endif()
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/plugins/roughconductor.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/plugins/roughconductor.so")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/plugins/roughconductor.so"
         RPATH "\$ORIGIN:\$ORIGIN/../drjit")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin/plugins" TYPE SHARED_LIBRARY FILES "/home/matt/rayference/eradiate/ext/mitsuba/out/build/eradiate/plugins/roughconductor.so")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/plugins/roughconductor.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/plugins/roughconductor.so")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/plugins/roughconductor.so"
         OLD_RPATH "/home/matt/rayference/eradiate/ext/mitsuba/out/build/eradiate:"
         NEW_RPATH "\$ORIGIN:\$ORIGIN/../drjit")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/llvm-strip-11" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/plugins/roughconductor.so")
    endif()
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/plugins/roughdielectric.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/plugins/roughdielectric.so")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/plugins/roughdielectric.so"
         RPATH "\$ORIGIN:\$ORIGIN/../drjit")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin/plugins" TYPE SHARED_LIBRARY FILES "/home/matt/rayference/eradiate/ext/mitsuba/out/build/eradiate/plugins/roughdielectric.so")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/plugins/roughdielectric.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/plugins/roughdielectric.so")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/plugins/roughdielectric.so"
         OLD_RPATH "/home/matt/rayference/eradiate/ext/mitsuba/out/build/eradiate:"
         NEW_RPATH "\$ORIGIN:\$ORIGIN/../drjit")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/llvm-strip-11" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/plugins/roughdielectric.so")
    endif()
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/plugins/roughplastic.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/plugins/roughplastic.so")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/plugins/roughplastic.so"
         RPATH "\$ORIGIN:\$ORIGIN/../drjit")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin/plugins" TYPE SHARED_LIBRARY FILES "/home/matt/rayference/eradiate/ext/mitsuba/out/build/eradiate/plugins/roughplastic.so")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/plugins/roughplastic.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/plugins/roughplastic.so")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/plugins/roughplastic.so"
         OLD_RPATH "/home/matt/rayference/eradiate/ext/mitsuba/out/build/eradiate:"
         NEW_RPATH "\$ORIGIN:\$ORIGIN/../drjit")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/llvm-strip-11" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/plugins/roughplastic.so")
    endif()
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/plugins/thindielectric.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/plugins/thindielectric.so")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/plugins/thindielectric.so"
         RPATH "\$ORIGIN:\$ORIGIN/../drjit")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin/plugins" TYPE SHARED_LIBRARY FILES "/home/matt/rayference/eradiate/ext/mitsuba/out/build/eradiate/plugins/thindielectric.so")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/plugins/thindielectric.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/plugins/thindielectric.so")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/plugins/thindielectric.so"
         OLD_RPATH "/home/matt/rayference/eradiate/ext/mitsuba/out/build/eradiate:"
         NEW_RPATH "\$ORIGIN:\$ORIGIN/../drjit")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/llvm-strip-11" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/plugins/thindielectric.so")
    endif()
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/plugins/twosided.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/plugins/twosided.so")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/plugins/twosided.so"
         RPATH "\$ORIGIN:\$ORIGIN/../drjit")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin/plugins" TYPE SHARED_LIBRARY FILES "/home/matt/rayference/eradiate/ext/mitsuba/out/build/eradiate/plugins/twosided.so")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/plugins/twosided.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/plugins/twosided.so")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/plugins/twosided.so"
         OLD_RPATH "/home/matt/rayference/eradiate/ext/mitsuba/out/build/eradiate:"
         NEW_RPATH "\$ORIGIN:\$ORIGIN/../drjit")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/llvm-strip-11" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/plugins/twosided.so")
    endif()
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/plugins/polarizer.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/plugins/polarizer.so")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/plugins/polarizer.so"
         RPATH "\$ORIGIN:\$ORIGIN/../drjit")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin/plugins" TYPE SHARED_LIBRARY FILES "/home/matt/rayference/eradiate/ext/mitsuba/out/build/eradiate/plugins/polarizer.so")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/plugins/polarizer.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/plugins/polarizer.so")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/plugins/polarizer.so"
         OLD_RPATH "/home/matt/rayference/eradiate/ext/mitsuba/out/build/eradiate:"
         NEW_RPATH "\$ORIGIN:\$ORIGIN/../drjit")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/llvm-strip-11" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/plugins/polarizer.so")
    endif()
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/plugins/retarder.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/plugins/retarder.so")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/plugins/retarder.so"
         RPATH "\$ORIGIN:\$ORIGIN/../drjit")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin/plugins" TYPE SHARED_LIBRARY FILES "/home/matt/rayference/eradiate/ext/mitsuba/out/build/eradiate/plugins/retarder.so")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/plugins/retarder.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/plugins/retarder.so")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/plugins/retarder.so"
         OLD_RPATH "/home/matt/rayference/eradiate/ext/mitsuba/out/build/eradiate:"
         NEW_RPATH "\$ORIGIN:\$ORIGIN/../drjit")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/llvm-strip-11" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/plugins/retarder.so")
    endif()
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/plugins/circular.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/plugins/circular.so")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/plugins/circular.so"
         RPATH "\$ORIGIN:\$ORIGIN/../drjit")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin/plugins" TYPE SHARED_LIBRARY FILES "/home/matt/rayference/eradiate/ext/mitsuba/out/build/eradiate/plugins/circular.so")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/plugins/circular.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/plugins/circular.so")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/plugins/circular.so"
         OLD_RPATH "/home/matt/rayference/eradiate/ext/mitsuba/out/build/eradiate:"
         NEW_RPATH "\$ORIGIN:\$ORIGIN/../drjit")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/llvm-strip-11" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/plugins/circular.so")
    endif()
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/plugins/measured_polarized.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/plugins/measured_polarized.so")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/plugins/measured_polarized.so"
         RPATH "\$ORIGIN:\$ORIGIN/../drjit")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin/plugins" TYPE SHARED_LIBRARY FILES "/home/matt/rayference/eradiate/ext/mitsuba/out/build/eradiate/plugins/measured_polarized.so")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/plugins/measured_polarized.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/plugins/measured_polarized.so")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/plugins/measured_polarized.so"
         OLD_RPATH "/home/matt/rayference/eradiate/ext/mitsuba/out/build/eradiate:"
         NEW_RPATH "\$ORIGIN:\$ORIGIN/../drjit")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/llvm-strip-11" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/plugins/measured_polarized.so")
    endif()
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/plugins/pplastic.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/plugins/pplastic.so")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/plugins/pplastic.so"
         RPATH "\$ORIGIN:\$ORIGIN/../drjit")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin/plugins" TYPE SHARED_LIBRARY FILES "/home/matt/rayference/eradiate/ext/mitsuba/out/build/eradiate/plugins/pplastic.so")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/plugins/pplastic.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/plugins/pplastic.so")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/plugins/pplastic.so"
         OLD_RPATH "/home/matt/rayference/eradiate/ext/mitsuba/out/build/eradiate:"
         NEW_RPATH "\$ORIGIN:\$ORIGIN/../drjit")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/llvm-strip-11" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/plugins/pplastic.so")
    endif()
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/plugins/principled.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/plugins/principled.so")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/plugins/principled.so"
         RPATH "\$ORIGIN:\$ORIGIN/../drjit")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin/plugins" TYPE SHARED_LIBRARY FILES "/home/matt/rayference/eradiate/ext/mitsuba/out/build/eradiate/plugins/principled.so")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/plugins/principled.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/plugins/principled.so")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/plugins/principled.so"
         OLD_RPATH "/home/matt/rayference/eradiate/ext/mitsuba/out/build/eradiate:"
         NEW_RPATH "\$ORIGIN:\$ORIGIN/../drjit")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/llvm-strip-11" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/plugins/principled.so")
    endif()
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/plugins/principledthin.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/plugins/principledthin.so")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/plugins/principledthin.so"
         RPATH "\$ORIGIN:\$ORIGIN/../drjit")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin/plugins" TYPE SHARED_LIBRARY FILES "/home/matt/rayference/eradiate/ext/mitsuba/out/build/eradiate/plugins/principledthin.so")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/plugins/principledthin.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/plugins/principledthin.so")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/plugins/principledthin.so"
         OLD_RPATH "/home/matt/rayference/eradiate/ext/mitsuba/out/build/eradiate:"
         NEW_RPATH "\$ORIGIN:\$ORIGIN/../drjit")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/llvm-strip-11" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/plugins/principledthin.so")
    endif()
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
endif()

