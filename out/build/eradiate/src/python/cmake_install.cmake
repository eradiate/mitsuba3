# Install script for directory: /home/matt/rayference/eradiate/ext/mitsuba/src/python

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
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/mitsuba/scalar_mono.cpython-312-x86_64-linux-gnu.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/mitsuba/scalar_mono.cpython-312-x86_64-linux-gnu.so")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/mitsuba/scalar_mono.cpython-312-x86_64-linux-gnu.so"
         RPATH "\$ORIGIN:\$ORIGIN/../drjit")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/mitsuba" TYPE MODULE FILES "/home/matt/rayference/eradiate/ext/mitsuba/out/build/eradiate/python/mitsuba/scalar_mono.cpython-312-x86_64-linux-gnu.so")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/mitsuba/scalar_mono.cpython-312-x86_64-linux-gnu.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/mitsuba/scalar_mono.cpython-312-x86_64-linux-gnu.so")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/mitsuba/scalar_mono.cpython-312-x86_64-linux-gnu.so"
         OLD_RPATH "/home/matt/rayference/eradiate/ext/mitsuba/out/build/eradiate/python/mitsuba:/home/matt/rayference/eradiate/ext/mitsuba/out/build/eradiate:"
         NEW_RPATH "\$ORIGIN:\$ORIGIN/../drjit")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/llvm-strip-11" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/mitsuba/scalar_mono.cpython-312-x86_64-linux-gnu.so")
    endif()
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/mitsuba/scalar_mono_double.cpython-312-x86_64-linux-gnu.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/mitsuba/scalar_mono_double.cpython-312-x86_64-linux-gnu.so")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/mitsuba/scalar_mono_double.cpython-312-x86_64-linux-gnu.so"
         RPATH "\$ORIGIN:\$ORIGIN/../drjit")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/mitsuba" TYPE MODULE FILES "/home/matt/rayference/eradiate/ext/mitsuba/out/build/eradiate/python/mitsuba/scalar_mono_double.cpython-312-x86_64-linux-gnu.so")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/mitsuba/scalar_mono_double.cpython-312-x86_64-linux-gnu.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/mitsuba/scalar_mono_double.cpython-312-x86_64-linux-gnu.so")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/mitsuba/scalar_mono_double.cpython-312-x86_64-linux-gnu.so"
         OLD_RPATH "/home/matt/rayference/eradiate/ext/mitsuba/out/build/eradiate/python/mitsuba:/home/matt/rayference/eradiate/ext/mitsuba/out/build/eradiate:"
         NEW_RPATH "\$ORIGIN:\$ORIGIN/../drjit")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/llvm-strip-11" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/mitsuba/scalar_mono_double.cpython-312-x86_64-linux-gnu.so")
    endif()
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/mitsuba/scalar_mono_polarized.cpython-312-x86_64-linux-gnu.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/mitsuba/scalar_mono_polarized.cpython-312-x86_64-linux-gnu.so")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/mitsuba/scalar_mono_polarized.cpython-312-x86_64-linux-gnu.so"
         RPATH "\$ORIGIN:\$ORIGIN/../drjit")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/mitsuba" TYPE MODULE FILES "/home/matt/rayference/eradiate/ext/mitsuba/out/build/eradiate/python/mitsuba/scalar_mono_polarized.cpython-312-x86_64-linux-gnu.so")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/mitsuba/scalar_mono_polarized.cpython-312-x86_64-linux-gnu.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/mitsuba/scalar_mono_polarized.cpython-312-x86_64-linux-gnu.so")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/mitsuba/scalar_mono_polarized.cpython-312-x86_64-linux-gnu.so"
         OLD_RPATH "/home/matt/rayference/eradiate/ext/mitsuba/out/build/eradiate/python/mitsuba:/home/matt/rayference/eradiate/ext/mitsuba/out/build/eradiate:"
         NEW_RPATH "\$ORIGIN:\$ORIGIN/../drjit")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/llvm-strip-11" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/mitsuba/scalar_mono_polarized.cpython-312-x86_64-linux-gnu.so")
    endif()
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/mitsuba/scalar_mono_polarized_double.cpython-312-x86_64-linux-gnu.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/mitsuba/scalar_mono_polarized_double.cpython-312-x86_64-linux-gnu.so")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/mitsuba/scalar_mono_polarized_double.cpython-312-x86_64-linux-gnu.so"
         RPATH "\$ORIGIN:\$ORIGIN/../drjit")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/mitsuba" TYPE MODULE FILES "/home/matt/rayference/eradiate/ext/mitsuba/out/build/eradiate/python/mitsuba/scalar_mono_polarized_double.cpython-312-x86_64-linux-gnu.so")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/mitsuba/scalar_mono_polarized_double.cpython-312-x86_64-linux-gnu.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/mitsuba/scalar_mono_polarized_double.cpython-312-x86_64-linux-gnu.so")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/mitsuba/scalar_mono_polarized_double.cpython-312-x86_64-linux-gnu.so"
         OLD_RPATH "/home/matt/rayference/eradiate/ext/mitsuba/out/build/eradiate/python/mitsuba:/home/matt/rayference/eradiate/ext/mitsuba/out/build/eradiate:"
         NEW_RPATH "\$ORIGIN:\$ORIGIN/../drjit")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/llvm-strip-11" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/mitsuba/scalar_mono_polarized_double.cpython-312-x86_64-linux-gnu.so")
    endif()
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/mitsuba/scalar_rgb.cpython-312-x86_64-linux-gnu.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/mitsuba/scalar_rgb.cpython-312-x86_64-linux-gnu.so")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/mitsuba/scalar_rgb.cpython-312-x86_64-linux-gnu.so"
         RPATH "\$ORIGIN:\$ORIGIN/../drjit")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/mitsuba" TYPE MODULE FILES "/home/matt/rayference/eradiate/ext/mitsuba/out/build/eradiate/python/mitsuba/scalar_rgb.cpython-312-x86_64-linux-gnu.so")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/mitsuba/scalar_rgb.cpython-312-x86_64-linux-gnu.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/mitsuba/scalar_rgb.cpython-312-x86_64-linux-gnu.so")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/mitsuba/scalar_rgb.cpython-312-x86_64-linux-gnu.so"
         OLD_RPATH "/home/matt/rayference/eradiate/ext/mitsuba/out/build/eradiate/python/mitsuba:/home/matt/rayference/eradiate/ext/mitsuba/out/build/eradiate:"
         NEW_RPATH "\$ORIGIN:\$ORIGIN/../drjit")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/llvm-strip-11" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/mitsuba/scalar_rgb.cpython-312-x86_64-linux-gnu.so")
    endif()
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/mitsuba/scalar_spectral.cpython-312-x86_64-linux-gnu.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/mitsuba/scalar_spectral.cpython-312-x86_64-linux-gnu.so")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/mitsuba/scalar_spectral.cpython-312-x86_64-linux-gnu.so"
         RPATH "\$ORIGIN:\$ORIGIN/../drjit")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/mitsuba" TYPE MODULE FILES "/home/matt/rayference/eradiate/ext/mitsuba/out/build/eradiate/python/mitsuba/scalar_spectral.cpython-312-x86_64-linux-gnu.so")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/mitsuba/scalar_spectral.cpython-312-x86_64-linux-gnu.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/mitsuba/scalar_spectral.cpython-312-x86_64-linux-gnu.so")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/mitsuba/scalar_spectral.cpython-312-x86_64-linux-gnu.so"
         OLD_RPATH "/home/matt/rayference/eradiate/ext/mitsuba/out/build/eradiate/python/mitsuba:/home/matt/rayference/eradiate/ext/mitsuba/out/build/eradiate:"
         NEW_RPATH "\$ORIGIN:\$ORIGIN/../drjit")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/llvm-strip-11" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/mitsuba/scalar_spectral.cpython-312-x86_64-linux-gnu.so")
    endif()
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/mitsuba/llvm_ad_rgb.cpython-312-x86_64-linux-gnu.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/mitsuba/llvm_ad_rgb.cpython-312-x86_64-linux-gnu.so")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/mitsuba/llvm_ad_rgb.cpython-312-x86_64-linux-gnu.so"
         RPATH "\$ORIGIN:\$ORIGIN/../drjit")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/mitsuba" TYPE MODULE FILES "/home/matt/rayference/eradiate/ext/mitsuba/out/build/eradiate/python/mitsuba/llvm_ad_rgb.cpython-312-x86_64-linux-gnu.so")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/mitsuba/llvm_ad_rgb.cpython-312-x86_64-linux-gnu.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/mitsuba/llvm_ad_rgb.cpython-312-x86_64-linux-gnu.so")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/mitsuba/llvm_ad_rgb.cpython-312-x86_64-linux-gnu.so"
         OLD_RPATH "/home/matt/rayference/eradiate/ext/mitsuba/out/build/eradiate/python/mitsuba:/home/matt/rayference/eradiate/ext/mitsuba/out/build/eradiate:"
         NEW_RPATH "\$ORIGIN:\$ORIGIN/../drjit")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/llvm-strip-11" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/mitsuba/llvm_ad_rgb.cpython-312-x86_64-linux-gnu.so")
    endif()
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/mitsuba/llvm_ad_spectral_polarized.cpython-312-x86_64-linux-gnu.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/mitsuba/llvm_ad_spectral_polarized.cpython-312-x86_64-linux-gnu.so")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/mitsuba/llvm_ad_spectral_polarized.cpython-312-x86_64-linux-gnu.so"
         RPATH "\$ORIGIN:\$ORIGIN/../drjit")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/mitsuba" TYPE MODULE FILES "/home/matt/rayference/eradiate/ext/mitsuba/out/build/eradiate/python/mitsuba/llvm_ad_spectral_polarized.cpython-312-x86_64-linux-gnu.so")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/mitsuba/llvm_ad_spectral_polarized.cpython-312-x86_64-linux-gnu.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/mitsuba/llvm_ad_spectral_polarized.cpython-312-x86_64-linux-gnu.so")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/mitsuba/llvm_ad_spectral_polarized.cpython-312-x86_64-linux-gnu.so"
         OLD_RPATH "/home/matt/rayference/eradiate/ext/mitsuba/out/build/eradiate/python/mitsuba:/home/matt/rayference/eradiate/ext/mitsuba/out/build/eradiate:"
         NEW_RPATH "\$ORIGIN:\$ORIGIN/../drjit")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/llvm-strip-11" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/mitsuba/llvm_ad_spectral_polarized.cpython-312-x86_64-linux-gnu.so")
    endif()
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/mitsuba/libnanobind-drjit.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/mitsuba/libnanobind-drjit.so")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/mitsuba/libnanobind-drjit.so"
         RPATH "\$ORIGIN:\$ORIGIN/../drjit")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/mitsuba" TYPE SHARED_LIBRARY FILES "/home/matt/rayference/eradiate/ext/mitsuba/out/build/eradiate/python/mitsuba/libnanobind-drjit.so")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/mitsuba/libnanobind-drjit.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/mitsuba/libnanobind-drjit.so")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/mitsuba/libnanobind-drjit.so"
         OLD_RPATH "::::::::::::::::::::::::"
         NEW_RPATH "\$ORIGIN:\$ORIGIN/../drjit")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/llvm-strip-11" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/mitsuba/libnanobind-drjit.so")
    endif()
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/mitsuba/mitsuba_ext.cpython-312-x86_64-linux-gnu.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/mitsuba/mitsuba_ext.cpython-312-x86_64-linux-gnu.so")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/mitsuba/mitsuba_ext.cpython-312-x86_64-linux-gnu.so"
         RPATH "\$ORIGIN:\$ORIGIN/../drjit")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/mitsuba" TYPE MODULE FILES "/home/matt/rayference/eradiate/ext/mitsuba/out/build/eradiate/python/mitsuba/mitsuba_ext.cpython-312-x86_64-linux-gnu.so")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/mitsuba/mitsuba_ext.cpython-312-x86_64-linux-gnu.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/mitsuba/mitsuba_ext.cpython-312-x86_64-linux-gnu.so")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/mitsuba/mitsuba_ext.cpython-312-x86_64-linux-gnu.so"
         OLD_RPATH "/home/matt/rayference/eradiate/ext/mitsuba/out/build/eradiate/python/mitsuba:/home/matt/rayference/eradiate/ext/mitsuba/out/build/eradiate:"
         NEW_RPATH "\$ORIGIN:\$ORIGIN/../drjit")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/llvm-strip-11" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/mitsuba/mitsuba_ext.cpython-312-x86_64-linux-gnu.so")
    endif()
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/mitsuba/mitsuba_alias.cpython-312-x86_64-linux-gnu.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/mitsuba/mitsuba_alias.cpython-312-x86_64-linux-gnu.so")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/mitsuba/mitsuba_alias.cpython-312-x86_64-linux-gnu.so"
         RPATH "\$ORIGIN:\$ORIGIN/../drjit")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/mitsuba" TYPE MODULE FILES "/home/matt/rayference/eradiate/ext/mitsuba/out/build/eradiate/python/mitsuba/mitsuba_alias.cpython-312-x86_64-linux-gnu.so")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/mitsuba/mitsuba_alias.cpython-312-x86_64-linux-gnu.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/mitsuba/mitsuba_alias.cpython-312-x86_64-linux-gnu.so")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/mitsuba/mitsuba_alias.cpython-312-x86_64-linux-gnu.so"
         OLD_RPATH "/home/matt/rayference/eradiate/ext/mitsuba/out/build/eradiate/python/mitsuba:"
         NEW_RPATH "\$ORIGIN:\$ORIGIN/../drjit")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/llvm-strip-11" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/mitsuba/mitsuba_alias.cpython-312-x86_64-linux-gnu.so")
    endif()
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/mitsuba" TYPE FILE FILES "/home/matt/rayference/eradiate/ext/mitsuba/src/python/__init__.py")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/mitsuba/mitsuba_stubs" TYPE FILE FILES "/home/matt/rayference/eradiate/ext/mitsuba/src/python/mitsuba_stubs/__init__.py")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/mitsuba/python" TYPE FILE FILES "/home/matt/rayference/eradiate/ext/mitsuba/src/python/python/__init__.py")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/mitsuba/python/ad" TYPE FILE FILES "/home/matt/rayference/eradiate/ext/mitsuba/src/python/python/ad/__init__.py")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/mitsuba/python/ad" TYPE FILE FILES "/home/matt/rayference/eradiate/ext/mitsuba/src/python/python/ad/guiding.py")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/mitsuba/python/ad/integrators" TYPE FILE FILES "/home/matt/rayference/eradiate/ext/mitsuba/src/python/python/ad/integrators/__init__.py")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/mitsuba/python/ad/integrators" TYPE FILE FILES "/home/matt/rayference/eradiate/ext/mitsuba/src/python/python/ad/integrators/common.py")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/mitsuba/python/ad/integrators" TYPE FILE FILES "/home/matt/rayference/eradiate/ext/mitsuba/src/python/python/ad/integrators/direct_projective.py")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/mitsuba/python/ad/integrators" TYPE FILE FILES "/home/matt/rayference/eradiate/ext/mitsuba/src/python/python/ad/integrators/prb.py")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/mitsuba/python/ad/integrators" TYPE FILE FILES "/home/matt/rayference/eradiate/ext/mitsuba/src/python/python/ad/integrators/prb_basic.py")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/mitsuba/python/ad/integrators" TYPE FILE FILES "/home/matt/rayference/eradiate/ext/mitsuba/src/python/python/ad/integrators/prb_projective.py")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/mitsuba/python/ad/integrators" TYPE FILE FILES "/home/matt/rayference/eradiate/ext/mitsuba/src/python/python/ad/integrators/prbvolpath.py")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/mitsuba/python/ad" TYPE FILE FILES "/home/matt/rayference/eradiate/ext/mitsuba/src/python/python/ad/largesteps.py")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/mitsuba/python/ad" TYPE FILE FILES "/home/matt/rayference/eradiate/ext/mitsuba/src/python/python/ad/optimizers.py")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/mitsuba/python/ad" TYPE FILE FILES "/home/matt/rayference/eradiate/ext/mitsuba/src/python/python/ad/projective.py")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/mitsuba/python" TYPE FILE FILES "/home/matt/rayference/eradiate/ext/mitsuba/src/python/python/chi2.py")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/mitsuba/python" TYPE FILE FILES "/home/matt/rayference/eradiate/ext/mitsuba/src/python/python/cli.py")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/mitsuba/python" TYPE FILE FILES "/home/matt/rayference/eradiate/ext/mitsuba/src/python/python/math_py.py")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/mitsuba/python" TYPE FILE FILES "/home/matt/rayference/eradiate/ext/mitsuba/src/python/python/polvis.py")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/mitsuba/python" TYPE FILE FILES "/home/matt/rayference/eradiate/ext/mitsuba/src/python/python/sys_info.py")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/mitsuba/python/test" TYPE FILE FILES "/home/matt/rayference/eradiate/ext/mitsuba/src/python/python/test/__init__.py")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/mitsuba/python/test" TYPE FILE FILES "/home/matt/rayference/eradiate/ext/mitsuba/src/python/python/test/util.py")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/mitsuba/python/tests" TYPE FILE FILES "/home/matt/rayference/eradiate/ext/mitsuba/src/python/python/tests/test_largesteps.py")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/mitsuba/python/tests" TYPE FILE FILES "/home/matt/rayference/eradiate/ext/mitsuba/src/python/python/tests/test_util.py")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/mitsuba/python" TYPE FILE FILES "/home/matt/rayference/eradiate/ext/mitsuba/src/python/python/tonemap.py")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/mitsuba/python" TYPE FILE FILES "/home/matt/rayference/eradiate/ext/mitsuba/src/python/python/util.py")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/mitsuba/python" TYPE FILE FILES "/home/matt/rayference/eradiate/ext/mitsuba/src/python/python/xml.py")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/mitsuba" TYPE FILE FILES "/home/matt/rayference/eradiate/ext/mitsuba/out/build/eradiate/src/python/../../python/mitsuba/config.py")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/mitsuba" TYPE FILE FILES
    "/home/matt/rayference/eradiate/ext/mitsuba/out/build/eradiate/python/mitsuba/__init__.pyi"
    "/home/matt/rayference/eradiate/ext/mitsuba/out/build/eradiate/python/mitsuba/py.typed"
    "/home/matt/rayference/eradiate/ext/mitsuba/out/build/eradiate/python/mitsuba/math.pyi"
    "/home/matt/rayference/eradiate/ext/mitsuba/out/build/eradiate/python/mitsuba/quad.pyi"
    "/home/matt/rayference/eradiate/ext/mitsuba/out/build/eradiate/python/mitsuba/spline.pyi"
    "/home/matt/rayference/eradiate/ext/mitsuba/out/build/eradiate/python/mitsuba/mueller.pyi"
    "/home/matt/rayference/eradiate/ext/mitsuba/out/build/eradiate/python/mitsuba/warp.pyi"
    "/home/matt/rayference/eradiate/ext/mitsuba/out/build/eradiate/python/mitsuba/misc.pyi"
    "/home/matt/rayference/eradiate/ext/mitsuba/out/build/eradiate/python/mitsuba/detail.pyi"
    "/home/matt/rayference/eradiate/ext/mitsuba/out/build/eradiate/python/mitsuba/filesystem.pyi"
    "/home/matt/rayference/eradiate/ext/mitsuba/out/build/eradiate/python/mitsuba/python/ad.pyi"
    "/home/matt/rayference/eradiate/ext/mitsuba/out/build/eradiate/python/mitsuba/python/chi2.pyi"
    "/home/matt/rayference/eradiate/ext/mitsuba/out/build/eradiate/python/mitsuba/python/math_py.pyi"
    "/home/matt/rayference/eradiate/ext/mitsuba/out/build/eradiate/python/mitsuba/python/util.pyi"
    "/home/matt/rayference/eradiate/ext/mitsuba/out/build/eradiate/python/mitsuba/python/xml.pyi"
    )
endif()

