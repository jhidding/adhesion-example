# adhesion-example
This is a short astrophysical program showing howto use CGAL (www.cgal.org) to compute the adhesion model, describing the large-scale structure of the Universe.

### requirements
This program is written in C++, using the 2011 standard. In addition to a recent GCC, to compile and run this example you need to have the development versions of the folowing libraries installed:
  * CGAL
  * FFTW 3.0

### features
The features are kept minimal:
  * Gaussian random fields: The program creates its own GRF to test the adhesion code on. This is the part of the program that depends on FFTW.
  * Glass initial conditions: Initial conditions are interpolated to a glass. This glass is generated on the fly, using Lloyd iteration. The Lloyd iteration makes use of periodic Voronoi tessellations.
  * Computing adhesion via the regular triangulation.
  * Selecting parts of the resulting web structures by either sphere or half-plane.
  * Writing output in ASCII to either Stanford PLY format or Wavefront OBJ. In the last case the density of sheets is written as a texture coordinate. This can be used when loading the OBJ to Blender (www.blender.org).

### building
To build this program, there is a bash make script included. Just typing:
 
 > ./make all
 
in the project folder will compile all .cc files in ./src and put the .o in a separate directory ./obj. The program is then linked to the executable, which is put in the project folder. If any dependancies are located outside the default search path, you can edit the top lines of the make script to change CFLAGS and LDFLAGS accordingly.
