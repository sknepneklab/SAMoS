# NOTE: C++11 required!!!

While most of the code still uses the C++98 standard, in order to compile it a C++11 compiler is required.
Due to problems with the boost/functional/factory library when compiling on Mac OS X (OS X High Sierra and later)
SAMoS no longer uses that library. Instead a light-weight class factory was implemented using the
C++11 standard (it relies on variadic templates).

# Soft Active Matter on Surfaces (SAMoS)

## 1. OS

The code has been developed and tested under Linux and Mac OS X. 
In principle it should be possible to build and run the code on 
a Windows machine. However, this would require modifying the CMake
build scripts. This has not been tested. 

## 2. TUTORIALS

A detailed tutorial on how to set up a simulation in SAMoS can be found in 
/path/to/SAMoS_install/doc/tutorial/tutorial.html

here */path/to/SAMoS_install/* is the directory where SAMoS is installed (see section 5 below).

## 3. REQUIREMENTS

SAMoS code is written in C++. Data analysis and initial configuration 
building tools are written in Python 2.7 and require NumPy to run.

* Modern C++ compiler supporting the C++11 standard 
* Boost libraries (1.48 or newer, in particular Spirit parser)
* GNU Scientific Library (GSL) - version 1.13 or newer 
* CMake (2.8 or newer) - it is recommended to install ccmake GUI
* Doxygen - optional but recommended (LaTeX support for building the PDF reference manual
* VTK library (version 5 or 6)
* CGAL library (version 4.3 or newer). NOTE: Code will fail to compile with CGAL 4.2 or older

**NOTE:** On Mac OS X, it is suggested to use Mac Ports to install all necessary libraries.

**NOTE:** SAMoS is able to generate VTP files as its output. We suggest installing and using ParaView 
to visualise the results. 

## 4. COMPILING 

a) Clone the code from the GitHub repository using 

    git clone https://github.com/sknepneklab/SAMoS.git

b) cd SAMoS (or the name of the directory you chose to download the code into) <br />
c) mkdir build <br />
d) cd build <br />
e) ccmake ../ <br />
f) Use the CMake's GUI to chose appropriate settings <br />
g) Press 'c' key several times to configure <br />
h) If all libraries have been found, CMake will allow you to create a Makefile. If not, please exit CMake GUI and install missing libraries. <br />
i) Press 'g' to create Make files. This will terminate CMake's GUI is return you to the shell. <br />
j) Type 'make -j 8' (-j option tells make how many parallel threads to use to compile; on a 4-core machine one can typically use 8 threads). <br /> 
k) If the compilation is successful, an executable 'samos' should appear in the build directory.

**NOTE:** Some Linux distributions with newer C++ compilers may have problems with running 8 parallel threads and can cause the machine to crash. 
If this happens, use '-j 2', which is likely not to cause any problems on most modern computers. 

**NOTE:** SAMoS uses many templated libraries in Boost. Compiling it may take several minutes even on a very fast machine.

**NOTE:** Depending on the compiler, you may get a number of warning messages. Those are harmless and you may safely ignore them.  

**NOTE:** SAMoS will only compile using c++11 (or newer) standard.

**WARNING:** Some Mac users have reported problems with compiling SAMoS on Sierra and High Sierra using packages installed with 'brew'. Switching to 
the 'MacPorts' seems to solve the problem.   

## 5. INSTALLING SAMoS 

After build has been successfully completed you can install SAMoS by typing 

make install 

in the build directory. This will install SAMoS package into $HOME/samos directory, where $HOME is 
the environment variable containing full path to your home directory. 

  - SAMoS binary (samos) will be placed in $HOME/samos/bin 
  - Examples will be in  $HOME/samos/examples 
  - Analysis scripts will be placed in $HOME/samos/analysis 
  - Basic tutorial files will be placed in $HOME/samos/doc/tutorial

Please make sure to add $HOME/samos/analysis to your PYTHONPATH shell variable and $HOME/samos/bin to the PATH variable.

**NOTE:** You can change the default installation directory by setting CMAKE_INSTALL_PREFIX variable in step 4. 


## 6. RUNNING SAMoS

The code requires two files to run:

   1. conf file containing simulation parameters, force field, integrator type, etc.
   2. data file containing initial position of particles (see configurations directory for a number of examples)

**NOTE:** Format of the configuration and data files can be found in the 'configurations' directory.

code is executed with 

./samos conf_file.conf

## 7. SOURCE DIRECTORY STRUCTURE

```
samos 
   /FormerAnalysis   - some earlier versions of scripts for data analysis
   /analysis         - current set of tools for analysing simulation results 
   /build            - build directory (contains the executable)
   /configurations   - contains set of directories with examples of different systems that can be studies with SAMoS. Some directories 
                       contain Python scripts for generating initial configurations. 
   /doc              - Doxygen files for generating documentation
   /utils            - Several additional utilities for building initial configuration 
   /src              - Source code 
      \aligner       - Implementations of different alignment interactions that control particle orientation 
         \pair       - Pairwise aligners (alignments between pairs of particles)
         \external   - Single particle aligners (such as alignment to the extrenal field)
      /constraints   - Implementations of constraints to different curved surfaces 
      /dump          - Handles output of the simulation data (it supports several standard formats, plain text, VTP, mol2, dcd, etc.)
      /integrators   - Implements several integrators of the equations of motion
      /log           - Handles logging of the simulation state (e.g., current time step, total energy, etc.)
      /messenger     - Handles info, warning and error messages produced by the code, as well as generating metadata (JSON of XML format)
      /parser        - Set of Boost Spirit parsers for parsing configuration files
      /population    - Handles 'population' control (e.g., cell division and death, particle type change, particle removal and addition, etc.)
      /potential     - Implements various interactions (some of them may not be actual potential)
         /angle      - Angle potentials for filament simulations
         /bond       - Bond potentials for filament simulations 
         /external   - Potentials (and forces) acting on a single particle (such as external field)
         /pair       - Pair (or multibody) forces acting as the result of interparticle interactions (e.g., Lennard-Jones)
      /system        - Definition of the base classes that define the system (particles, simulation box, mesh, etc.)
      /utils         - Several utility functions (e.g., random number generator classes)
```

## 8. CREDITS

### Lead developer:

Rastko Sknepnek (University of Dundee, UK) 

### Major contributors to the code:

Silke Henkes (University of Aberdeen, UK)   - many tools for building inital configurations and analysis tools <br />
Daniel Barton (University of Dundee, UK)    - tools for building and analysing tissue mechanics simulations <br />
Amit Das (National Institute for Biological Sciences, India)  - tools for building and analysing actomyosin simulations 




