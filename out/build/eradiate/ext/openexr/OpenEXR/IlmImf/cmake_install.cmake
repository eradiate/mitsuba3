# Install script for directory: /home/matt/rayference/eradiate/ext/mitsuba/ext/openexr/OpenEXR/IlmImf

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
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libIlmImf-mitsuba.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libIlmImf-mitsuba.so")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libIlmImf-mitsuba.so"
         RPATH "$ORIGIN:$ORIGIN/../drjit")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE SHARED_LIBRARY FILES "/home/matt/rayference/eradiate/ext/mitsuba/out/build/eradiate/libIlmImf-mitsuba.so")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libIlmImf-mitsuba.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libIlmImf-mitsuba.so")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libIlmImf-mitsuba.so"
         OLD_RPATH "/home/matt/rayference/eradiate/ext/mitsuba/out/build/eradiate:"
         NEW_RPATH "$ORIGIN:$ORIGIN/../drjit")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/llvm-strip-11" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libIlmImf-mitsuba.so")
    endif()
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/OpenEXR" TYPE FILE FILES
    "/home/matt/rayference/eradiate/ext/mitsuba/ext/openexr/OpenEXR/IlmImf/ImfForward.h"
    "/home/matt/rayference/eradiate/ext/mitsuba/ext/openexr/OpenEXR/IlmImf/ImfExport.h"
    "/home/matt/rayference/eradiate/ext/mitsuba/ext/openexr/OpenEXR/IlmImf/ImfAttribute.h"
    "/home/matt/rayference/eradiate/ext/mitsuba/ext/openexr/OpenEXR/IlmImf/ImfBoxAttribute.h"
    "/home/matt/rayference/eradiate/ext/mitsuba/ext/openexr/OpenEXR/IlmImf/ImfCRgbaFile.h"
    "/home/matt/rayference/eradiate/ext/mitsuba/ext/openexr/OpenEXR/IlmImf/ImfChannelList.h"
    "/home/matt/rayference/eradiate/ext/mitsuba/ext/openexr/OpenEXR/IlmImf/ImfChannelListAttribute.h"
    "/home/matt/rayference/eradiate/ext/mitsuba/ext/openexr/OpenEXR/IlmImf/ImfCompressionAttribute.h"
    "/home/matt/rayference/eradiate/ext/mitsuba/ext/openexr/OpenEXR/IlmImf/ImfDoubleAttribute.h"
    "/home/matt/rayference/eradiate/ext/mitsuba/ext/openexr/OpenEXR/IlmImf/ImfFloatAttribute.h"
    "/home/matt/rayference/eradiate/ext/mitsuba/ext/openexr/OpenEXR/IlmImf/ImfFrameBuffer.h"
    "/home/matt/rayference/eradiate/ext/mitsuba/ext/openexr/OpenEXR/IlmImf/ImfHeader.h"
    "/home/matt/rayference/eradiate/ext/mitsuba/ext/openexr/OpenEXR/IlmImf/ImfIO.h"
    "/home/matt/rayference/eradiate/ext/mitsuba/ext/openexr/OpenEXR/IlmImf/ImfInputFile.h"
    "/home/matt/rayference/eradiate/ext/mitsuba/ext/openexr/OpenEXR/IlmImf/ImfIntAttribute.h"
    "/home/matt/rayference/eradiate/ext/mitsuba/ext/openexr/OpenEXR/IlmImf/ImfLineOrderAttribute.h"
    "/home/matt/rayference/eradiate/ext/mitsuba/ext/openexr/OpenEXR/IlmImf/ImfMatrixAttribute.h"
    "/home/matt/rayference/eradiate/ext/mitsuba/ext/openexr/OpenEXR/IlmImf/ImfOpaqueAttribute.h"
    "/home/matt/rayference/eradiate/ext/mitsuba/ext/openexr/OpenEXR/IlmImf/ImfOutputFile.h"
    "/home/matt/rayference/eradiate/ext/mitsuba/ext/openexr/OpenEXR/IlmImf/ImfRgbaFile.h"
    "/home/matt/rayference/eradiate/ext/mitsuba/ext/openexr/OpenEXR/IlmImf/ImfStringAttribute.h"
    "/home/matt/rayference/eradiate/ext/mitsuba/ext/openexr/OpenEXR/IlmImf/ImfVecAttribute.h"
    "/home/matt/rayference/eradiate/ext/mitsuba/ext/openexr/OpenEXR/IlmImf/ImfHuf.h"
    "/home/matt/rayference/eradiate/ext/mitsuba/ext/openexr/OpenEXR/IlmImf/ImfWav.h"
    "/home/matt/rayference/eradiate/ext/mitsuba/ext/openexr/OpenEXR/IlmImf/ImfLut.h"
    "/home/matt/rayference/eradiate/ext/mitsuba/ext/openexr/OpenEXR/IlmImf/ImfArray.h"
    "/home/matt/rayference/eradiate/ext/mitsuba/ext/openexr/OpenEXR/IlmImf/ImfCompression.h"
    "/home/matt/rayference/eradiate/ext/mitsuba/ext/openexr/OpenEXR/IlmImf/ImfLineOrder.h"
    "/home/matt/rayference/eradiate/ext/mitsuba/ext/openexr/OpenEXR/IlmImf/ImfName.h"
    "/home/matt/rayference/eradiate/ext/mitsuba/ext/openexr/OpenEXR/IlmImf/ImfPixelType.h"
    "/home/matt/rayference/eradiate/ext/mitsuba/ext/openexr/OpenEXR/IlmImf/ImfVersion.h"
    "/home/matt/rayference/eradiate/ext/mitsuba/ext/openexr/OpenEXR/IlmImf/ImfXdr.h"
    "/home/matt/rayference/eradiate/ext/mitsuba/ext/openexr/OpenEXR/IlmImf/ImfConvert.h"
    "/home/matt/rayference/eradiate/ext/mitsuba/ext/openexr/OpenEXR/IlmImf/ImfPreviewImage.h"
    "/home/matt/rayference/eradiate/ext/mitsuba/ext/openexr/OpenEXR/IlmImf/ImfPreviewImageAttribute.h"
    "/home/matt/rayference/eradiate/ext/mitsuba/ext/openexr/OpenEXR/IlmImf/ImfChromaticities.h"
    "/home/matt/rayference/eradiate/ext/mitsuba/ext/openexr/OpenEXR/IlmImf/ImfChromaticitiesAttribute.h"
    "/home/matt/rayference/eradiate/ext/mitsuba/ext/openexr/OpenEXR/IlmImf/ImfKeyCode.h"
    "/home/matt/rayference/eradiate/ext/mitsuba/ext/openexr/OpenEXR/IlmImf/ImfKeyCodeAttribute.h"
    "/home/matt/rayference/eradiate/ext/mitsuba/ext/openexr/OpenEXR/IlmImf/ImfTimeCode.h"
    "/home/matt/rayference/eradiate/ext/mitsuba/ext/openexr/OpenEXR/IlmImf/ImfTimeCodeAttribute.h"
    "/home/matt/rayference/eradiate/ext/mitsuba/ext/openexr/OpenEXR/IlmImf/ImfRational.h"
    "/home/matt/rayference/eradiate/ext/mitsuba/ext/openexr/OpenEXR/IlmImf/ImfRationalAttribute.h"
    "/home/matt/rayference/eradiate/ext/mitsuba/ext/openexr/OpenEXR/IlmImf/ImfFramesPerSecond.h"
    "/home/matt/rayference/eradiate/ext/mitsuba/ext/openexr/OpenEXR/IlmImf/ImfStandardAttributes.h"
    "/home/matt/rayference/eradiate/ext/mitsuba/ext/openexr/OpenEXR/IlmImf/ImfStdIO.h"
    "/home/matt/rayference/eradiate/ext/mitsuba/ext/openexr/OpenEXR/IlmImf/ImfEnvmap.h"
    "/home/matt/rayference/eradiate/ext/mitsuba/ext/openexr/OpenEXR/IlmImf/ImfEnvmapAttribute.h"
    "/home/matt/rayference/eradiate/ext/mitsuba/ext/openexr/OpenEXR/IlmImf/ImfInt64.h"
    "/home/matt/rayference/eradiate/ext/mitsuba/ext/openexr/OpenEXR/IlmImf/ImfRgba.h"
    "/home/matt/rayference/eradiate/ext/mitsuba/ext/openexr/OpenEXR/IlmImf/ImfTileDescription.h"
    "/home/matt/rayference/eradiate/ext/mitsuba/ext/openexr/OpenEXR/IlmImf/ImfTileDescriptionAttribute.h"
    "/home/matt/rayference/eradiate/ext/mitsuba/ext/openexr/OpenEXR/IlmImf/ImfTiledInputFile.h"
    "/home/matt/rayference/eradiate/ext/mitsuba/ext/openexr/OpenEXR/IlmImf/ImfTiledOutputFile.h"
    "/home/matt/rayference/eradiate/ext/mitsuba/ext/openexr/OpenEXR/IlmImf/ImfTiledRgbaFile.h"
    "/home/matt/rayference/eradiate/ext/mitsuba/ext/openexr/OpenEXR/IlmImf/ImfRgbaYca.h"
    "/home/matt/rayference/eradiate/ext/mitsuba/ext/openexr/OpenEXR/IlmImf/ImfTestFile.h"
    "/home/matt/rayference/eradiate/ext/mitsuba/ext/openexr/OpenEXR/IlmImf/ImfThreading.h"
    "/home/matt/rayference/eradiate/ext/mitsuba/ext/openexr/OpenEXR/IlmImf/ImfB44Compressor.h"
    "/home/matt/rayference/eradiate/ext/mitsuba/ext/openexr/OpenEXR/IlmImf/ImfStringVectorAttribute.h"
    "/home/matt/rayference/eradiate/ext/mitsuba/ext/openexr/OpenEXR/IlmImf/ImfMultiView.h"
    "/home/matt/rayference/eradiate/ext/mitsuba/ext/openexr/OpenEXR/IlmImf/ImfAcesFile.h"
    "/home/matt/rayference/eradiate/ext/mitsuba/ext/openexr/OpenEXR/IlmImf/ImfMultiPartOutputFile.h"
    "/home/matt/rayference/eradiate/ext/mitsuba/ext/openexr/OpenEXR/IlmImf/ImfGenericOutputFile.h"
    "/home/matt/rayference/eradiate/ext/mitsuba/ext/openexr/OpenEXR/IlmImf/ImfMultiPartInputFile.h"
    "/home/matt/rayference/eradiate/ext/mitsuba/ext/openexr/OpenEXR/IlmImf/ImfGenericInputFile.h"
    "/home/matt/rayference/eradiate/ext/mitsuba/ext/openexr/OpenEXR/IlmImf/ImfPartType.h"
    "/home/matt/rayference/eradiate/ext/mitsuba/ext/openexr/OpenEXR/IlmImf/ImfPartHelper.h"
    "/home/matt/rayference/eradiate/ext/mitsuba/ext/openexr/OpenEXR/IlmImf/ImfOutputPart.h"
    "/home/matt/rayference/eradiate/ext/mitsuba/ext/openexr/OpenEXR/IlmImf/ImfTiledOutputPart.h"
    "/home/matt/rayference/eradiate/ext/mitsuba/ext/openexr/OpenEXR/IlmImf/ImfInputPart.h"
    "/home/matt/rayference/eradiate/ext/mitsuba/ext/openexr/OpenEXR/IlmImf/ImfTiledInputPart.h"
    "/home/matt/rayference/eradiate/ext/mitsuba/ext/openexr/OpenEXR/IlmImf/ImfDeepScanLineOutputFile.h"
    "/home/matt/rayference/eradiate/ext/mitsuba/ext/openexr/OpenEXR/IlmImf/ImfDeepScanLineOutputPart.h"
    "/home/matt/rayference/eradiate/ext/mitsuba/ext/openexr/OpenEXR/IlmImf/ImfDeepScanLineInputFile.h"
    "/home/matt/rayference/eradiate/ext/mitsuba/ext/openexr/OpenEXR/IlmImf/ImfDeepScanLineInputPart.h"
    "/home/matt/rayference/eradiate/ext/mitsuba/ext/openexr/OpenEXR/IlmImf/ImfDeepTiledInputFile.h"
    "/home/matt/rayference/eradiate/ext/mitsuba/ext/openexr/OpenEXR/IlmImf/ImfDeepTiledInputPart.h"
    "/home/matt/rayference/eradiate/ext/mitsuba/ext/openexr/OpenEXR/IlmImf/ImfDeepTiledOutputFile.h"
    "/home/matt/rayference/eradiate/ext/mitsuba/ext/openexr/OpenEXR/IlmImf/ImfDeepTiledOutputPart.h"
    "/home/matt/rayference/eradiate/ext/mitsuba/ext/openexr/OpenEXR/IlmImf/ImfDeepFrameBuffer.h"
    "/home/matt/rayference/eradiate/ext/mitsuba/ext/openexr/OpenEXR/IlmImf/ImfDeepCompositing.h"
    "/home/matt/rayference/eradiate/ext/mitsuba/ext/openexr/OpenEXR/IlmImf/ImfCompositeDeepScanLine.h"
    "/home/matt/rayference/eradiate/ext/mitsuba/ext/openexr/OpenEXR/IlmImf/ImfNamespace.h"
    "/home/matt/rayference/eradiate/ext/mitsuba/ext/openexr/OpenEXR/IlmImf/ImfDeepImageState.h"
    "/home/matt/rayference/eradiate/ext/mitsuba/ext/openexr/OpenEXR/IlmImf/ImfDeepImageStateAttribute.h"
    "/home/matt/rayference/eradiate/ext/mitsuba/ext/openexr/OpenEXR/IlmImf/ImfFloatVectorAttribute.h"
    )
endif()

