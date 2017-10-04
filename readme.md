# MuSTEM

MuSTEM is a transmission electron microscopy (TEM) simulation suite, in particular for scanning transmission electron microscopy (STEM) images, that was developed mainly at the University of Melbourne. The computing suite is based on the multislice
method.

Thermal scattering is accounted for with one of two models:

1. The absorptive scattering potential approach [1, 2, 7].
2. The quantum excitation of phonons (QEP) model [6] which provides overall results numerically equivalent to the frozen phonon (FPh) method [8] but provides different physical insights and allows both elastic and inelastic phonon scattering to be separately elucidated.

Image simulation of inner shell ionization, i.e. electron energy loss spectroscopy (EELS) and energy dispersive x-ray spectroscopy (EDX), are based on a parameterization of the effective ionization potential in the local approximation. For EELS, a correction is made for the finite detector aperture size. It is assumed that the aperture size is large enough that the elastic scattering of the outgoing electron can be neglected. More detail can be found in the following scientific paper:

[Modelling the inelastic scattering of fast electrons,
L.J. Allen, A.J. D'Alfonso and S.D. Findlay,
Ultramicroscopy, Vol. 151, pp. 11-22, (2015).](http://www.sciencedirect.com/science/article/pii/S0304399114002034)

## Getting Started

MuSTEM is built using the [PGI Fortran compiler and Microsoft Visual Studio 2015](https://www.pgroup.com/products/pvf.htm), the following arguements are passed to the compiler:

Build commands:

-Mpreprocess -Dsingle_precision -Bstatic -Mbackslash -mp -Mcuda=cuda8.0 -I"C:\Program Files\PGI\win64\17.3\include" -I"c:\program files\pgi\win64\17.3\include" -I"C:\Program Files\PGI\Microsoft Open Tools 14\include" -I"C:\Program Files (x86)\Windows Kits\10\Include\shared" -I"C:\Program Files (x86)\Windows Kits\10\Include\um" -fast -ta=tesla -Minform=warn 

Linker commands:

-Bstatic -mp -Mcuda=cuda8.0 -ta=tesla -o"$(Outdir)\MuSTEM_Open.exe" -Wl,/libpath:"C:\Program Files\PGI\win64\17.3\lib" -Wl,/libpath:"C:\Program Files\PGI\win64\2017\cuda\8.0\lib64" cufft.lib 

The links to the PGI CUDA libraries and Windows kits may need to be modified depending on the install directories. $(Outdir) represents the output directory of the Microsoft Visual Studio build.

### Prerequisites

[PGI Visual Fortran and Microsoft Visual Studio](https://www.pgroup.com/products/pvf.htm)

```
Give examples
```

### Installing

A step by step series of examples that tell you have to get a development env running

Say what the step will be

```
Give the example
```

And repeat

```
until finished
```

End with an example of getting some data out of the system or using it for a little demo

## Running the tests

Explain how to run the automated tests for this system

### Break down into end to end tests

Explain what these tests test and why

```
Give an example
```

### And coding style tests

Explain what these tests test and why

```
Give an example
```

## Deployment

Add additional notes about how to deploy this on a live system

## Built With

* [Dropwizard](http://www.dropwizard.io/1.0.2/docs/) - The web framework used
* [Maven](https://maven.apache.org/) - Dependency Management
* [ROME](https://rometools.github.io/rome/) - Used to generate RSS Feeds

## Contributing

Please read [CONTRIBUTING.md](https://gist.github.com/PurpleBooth/b24679402957c63ec426) for details on our code of conduct, and the process for submitting pull requests to us.

## Versioning

We use [SemVer](http://semver.org/) for versioning. For the versions available, see the [tags on this repository](https://github.com/your/project/tags). 

## Authors

* **Billie Thompson** - *Initial work* - [PurpleBooth](https://github.com/PurpleBooth)

## License

This project is licensed under the GNU GPLv3.0 - see the [LICENSE.md](LICENSE.md) file for details


## Acknowledgments

Les Allen, Adrian D'Alfonso and Scott Findlay originally took the initiative to make this code publicly available, with Adrian D'Alfonso taking on the considerable task of setting up the GPU code. Hamish Brown and Ben Forbes have subsequently made substantial refinements and extensions, with Ben Forbes responsible for several efficiency improvements. The code was developed mainly at the University of Melbourne. We would like to acknowledge the contributions of Mark Oxley and Chris Rossouw to the theoretical and numerical aspects of collaborative research that underpins μSTEM 4.9.

In particular, the code IONIZER, whose primary author has been Mark Oxley, was used to set up the parametrized atomic scattering factors for ionization used in μSTEM 4.9. Eireann Cosgriff, Torgny Josefsson, Nathan Lugg, Andrew Martin, Gary Ruben and Chris Witte have been or still are members of the group at Melbourne University at various stages and all made contributions to our understanding of inelastic scattering and the simulation thereof. 

