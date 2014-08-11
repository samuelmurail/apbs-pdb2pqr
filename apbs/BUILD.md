APBS 1.4.1-binary BUILD:tada::clap:
=====================================

Contained herein are the steps that went in to creating the 1.4.1-binary
APBS release, and the branch upon which the release is based.

First, a few notes are in order.

- Only 64-bit operating systems are supported with this release.
- [FETK](http://www.fetk.org) support is included for the OS X and Linux builds, but not for the Windows build.  This is mostly to do with the inability of Windows and Autotools to play nicely together.  _NB: FETK has been modified somewhat to better inter-operate with APBS._
- Intel C and C++ compilers (version 14.0.3.20140415) and libraries were used for this release.
- CMake is required, as well as C and C++ compilers.
- If building the documentation, note that there is a bug in Doxygen that causes PDF generation to fail (see https://github.com/doxygen/doxygen/pull/178).


OS-Specific Steps
-----------------
Below are the specific steps used to produce this APBS release.  Following the steps below for your specific OS should produce the expected results.  However, a binary is not guaranteed: your mileage may vary.

In the following sections _&lt;INSTALL_DIR&gt;_ is the location to which the compiled and linked files will be copied.  <APBS> is the directory that was cloned from the Git repository (normally apbs-pdb2pqr).


OS X
----
The following steps were performed on OS X 10.9.4.  Building the executable is fairly straightforward:

1. Create a directory within which to build APBS.

    ````
    cd <APBS>/apbs
    mkdir osx
    cd osx
    ````
4. Build FETK.

    ```
    ditto ../fetk-1.5_APBS fetk
    cd fetk
    CC=icc CXX=icpc ./fetk-build maloc punc gamer mc
    cd ..
    ```

8. Run CMake to configure the build.

    ```
    cmake -DCMAKE_C_COMPILER=icc -DCMAKE_CXX_COMPILER=icpc -DCMAKE_BUILD_TYPE=RelWithDebInfo -DENABLE_FETK=ON \
    -DFETK_PATH=./fetk/build/x86_64-apple-darwin13.3.0 -DBUILD_SHARED_LIBS=OFF -DENABLE_BEM=OFF \
    -DBUILD_DOC=ON -DCMAKE_INSTALL_PREFIX=<INSTALL_DIR> ..
    ```
9. Compile and install!

    ```
    cmake --build . --target install
    ```

The result however isn't quite right.  The Intel compiler produces a binary with OpenMP support, and it relies on an Intel shared library.  However, the path to the shared library isn't baked into the binary.  Furthermore, the shared library paths that are baked into the binary are pointing to the wrong place for distribution.

The Intel compiler ships with, and links against Intel's OpenMP library libiomp5.  By default it's located in "/opt/intel/lib".  The APBS binary will not run unless it can find libiomp5.dylib, either by setting the DYLD_LIBRARY_PATH, or DYLD_FALLBACK_LIBRARY_PATH, or by "baking" the location into the binary.

Everything is now located in one place, but we still need to update the binary to point to the libraries that we built (and the missing Intel library) using a relative path.  Before we do that though, it would be nice to package this up in an app bundle, which is expected on the platform.

Eventually we'll utilize CMake to help with this task.  For now however, we'll mirror what was previously done for Windows and add some OS-specific supplementary material to `<APBS>/tools/osx` to help with the process.  This is still very much a manual process however.  (Yes a "trivial" shell script could take care of this, and yes it would imply a semblance of permanence which I still want to avoid.)


1. The plan is to create the App Bundle in the same directory that was used with CMake to build the installation image.  First we'll create the necessary directory structure.

    ```
    cd <INSTALL_DIR>
    mkdir -p APBS.app/Contents/MacOS
    mkdir -p APBS.app/Contents/Resources
    mkdir -p APBS.app/Contents/Frameworks
    ```

2. Next we need to copy the APBS binary and requisite libraries into the bundle.

    ```
    cd APBS.app/Contents
    cp <APBS>/tools/osx/Info.plist .
    cp <APBS>/tools/osx/apbs_term MacOS
    cp ../../bin/apbs MacOS
    ditto ../../share/apbs/tools/bin MacOS
    ditto <APBS>/apbs/osx/fetk/build/x86_64-apple-darwin13.3.0/lib Frameworks
    cp /opt/intel/lib/libiomp5.dylib Frameworks
    ```

3. Note that `otool -L MacOS/apbs` produces a list of dynamic library paths against which the binary is linked.  The binary needs to be updated to point to the libs in the "Frameworks" directory, as do all all of the other binaries in the MacOS directory.  This is accomplished using the _install_name_tool_ command.  A simple script to accomplish this has been created and interred in <APBS>/tools/osx.  Simply run the script, and the binary files will reference the dylibs in the Frameworks directory.

    ```
    <APBS>/tools/osx/patch_binary.sh
    ```

Linux
-----
Note that, at least on Ubuntu 14.04 LTS, in addition to cmake and the gnu or Intel compilers, you will also need to install libreadline-dev.

For linux, the build is the same as OS X:

1. Create a directory within which to build APBS.

    ````
    cd <APBS>/apbs
    mkdir linux
    cd linux
    ````
4. Build FETK.

    ```
    cp -r ../fetk-1.5_APBS fetk
    cd fetk
    LD=xild AR=xiar CC=icc CXX=icpc ./fetk-build maloc punc gamer mc
    cd ..
    ```

8. Run CMake to configure the build.

    ```
    cmake -DCMAKE_C_COMPILER=icc -DCMAKE_CXX_COMPILER=icpc -DCMAKE_BUILD_TYPE=RelWithDebInfo -DENABLE_FETK=ON \
    -DFETK_PATH=./fetk/build/x86_64-unknown-linux-gnu -DBUILD_SHARED_LIBS=OFF -DENABLE_BEM=OFF \
    -DBUILD_DOC=ON -DCMAKE_INSTALL_PREFIX=<INSTALL_DIR> ..
    ```
9. Compile and install!

    ```
    cmake --build . --target install
    ```

The results may then be installed in their usual locations.  Note that the directory structure created by CMake is conducive to this.

Windows 7
---------
Windows, requires additional "special" treatment.  To wit:

- Additional changes to CMakeLists.txt were required to Do the Right Thing with regard to Windows linker needs.
- Since Windows doesn't get the full FETK-love (due to no Win-autotools) we use the previously hacked upon MALOC that has been CMake-ified.  However, we want to use the Intel compiler, so we needed to forward those options.
- In two places in the source (src/generic/vhal.h and src/mg/vgrid.c) there were WIN32 #ifdef checks, however they needed to be _WIN32 to avoid great compiler-pain and ultimate (build) death.

After that, all that is necessary is to run CMake to configure the build:

```
cmake -G "NMake Makefiles" -DBUILD_SHARED_LIBS=OFF -DENABLE_BEM=OFF -DCMAKE_C_COMPILER=icl -DCMAKE_CXX_COMPILER=icl \
-DCMAKE_BUILD_TYPE=RelWithDebInfo -DCMAKE_INSTALL_PREFIX=<INSTALL_DIR> ..
```

And CMake to compile and install:

```
cmake --bulid . --target install
```

There exists an NSI configuration file in <APBS>/apbs/tools/windows/NSI Configuration that we use to create an installation program.
