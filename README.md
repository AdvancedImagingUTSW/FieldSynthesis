# FieldSynthesis

## Abstract

We introduce Field Synthesis, a theorem that can be used to synthesize any scanned or dithered light-sheet, including those used in lattice light-sheet microscopy (LLSM), from an incoherent superposition of one-dimensional intensity distributions. This user-friendly and modular approach offers a drastically simplified optical design, higher light-throughput, simultaneous multicolor illumination, and a 100% spatial duty cycle, thereby providing uncompromised biological imaging with decreased rates of photobleaching. 

## Manuscript


Bo-Jui Chang<sup>1</sup>, Mark Kittisopikul<sup>2,4</sup>, Kevin M. Dean<sup>1,3</sup>, Phillipe Roudot<sup>1,3</sup>, Erik Welf<sup>1,3</sup> and Reto Fiolka<sup>1,3</sup>.
"Democratizing Lattice Light-Sheet Microscopy with Field Synthesis."
 
### Affiliations
1. Department of Cell Biology, UT Southwestern Medical Center, Dallas, TX, USA. 
2. Department of Biophysics, UT Southwestern Medical Center, Dallas, TX, USA.
3. Lyda Hill Department of Bioinformatics, UT Southwestern Medical Center, Dallas, TX, USA.
4. Department of Cell and Molecular Biology, Feinberg School of Medicine, Northwestern University, Chicago, IL, USA.

## System Requirements

* [MATLAB (2017a or above), Mathworks, Natick, MA](https://www.mathworks.com/support/sysreq.html)
** Image Processing Toolbox
** Follow MATLAB link above for operating system requirements
* Git 2.18.0 or above

## Installation and Demo

Typical Install Time: 5 minutes

git clone https://github.com/AdvancedImagingUTSW/FieldSynthesis.git

See below documentation for demonstration. Typical run time: 10 minutes

## Code

### FieldSynthesisTheorem.m

   Small program to illustrate a new Field Synthesis Theorem.

   In essence it says that the projection of the absolute modulus of a
   complex field is the same as when one takes a sliding window in the
   Fourier domain, makes an inverse FFT of each slice, take the absolute
   modulus of that and sum it up while moving the window through the
   spectrum.  This has important applications for scanned light-sheets and
   how to generate them.

   Reto Fiolka, May 2017  
   Mark Kittisopikul, May 2017 - Aug 2018

   #### INPUT
   * efield - electric field at the focal plane, may be real or complex
            valued
   
   #### OUTPUT
   * efield - electric field at the focal plane
   * slice  - intensity of illumination pattern by field synthesis
   * smear  - intensity of illumination pattern by dithering
   * Q      - Fourier transform of individual line scan without phasing,a=10
   * T      - Fourier transform of individual line scan with phasing,a=10

```matlab
FieldSynthesisTheorem.m;
```

![Line scan simulation showing that field synthesis and lattice light sheet profiles are the same.](images/FieldSynthesisTheorem_1.png)
![Individual line scan analysis showing complex modulation depending on line scan position](images/FieldSynthesisTheorem_2.png)


### FieldSynthesisInteractive.m

FieldSynthesisInteractive Create an interactive line scan demonstration of
field synthesis

 #### INPUT
 * mask - mask at the pupil, which is the Fourier transform of electrical
        field at the focal plane
 * doshift - if true, shift the Fourier transform of the mask so the first
           pixel is in the center of the image rather than the upper left

 #### OUTPUT
 * hfig - handle for the display figure

 #### INTERACTIVE
 * The button in the lower left plays / pauses the movie.
 * The arrow buttons on the slider will move the scan by one column.
 * Clicking on the trough of the slider will move the scan by five columns.
 * The button in the lower right labeled R will reset the cumulative view.

 #### EXAMPLE
```matlab
 FieldSynthesisInteractive; % default demonstration with cameraman
```
![Screenshot of interactive field synthesis line scan simulation](images/FieldSynthesisInteractive.png)

```matlab
 FieldSynthesisInteractive(createAnnulus(),true); % demonstrate a Bessel beam 
```

![Screenshot of interactive field synthesis line scan simulation of Bessel beam](images/FieldSynthesisInteractive_bessel.png)

 Mark Kittisopikul , August 2018  
 Goldman Lab  
 Northwestern University

### FieldSynthesisVersusLattice.m

Simulation for field synthesis

   compares field synthesis vs square lattice


   Reto, May 2017  
   Mark Kittisopikul, August 2018

   #### INPUT
   * n - Defines the size of the image and mask to be n x n
   * w - Width of the mask components
   * r - Radius of the annulus (width is centered on the annulus)
   * offset - Offset of the side components of the square lattice
   * dispRange - Set which part of mask to display in figures

   #### OUTPUT
   * out - struct containing workspace of this function

```matlab
FieldSynthesisVersusLattice.m
```

![Profile demonstration showing identical profiles between field synthesis and dithered lattice light sheets](images/FieldSynthesisVersusLattice_Profiles.png)
![Demonstration showing Fourier components of lattice light sheet simulation](images/FieldSynthesisVersusLattice.png)


### createAnnulus.m

 #### INPUT (all optional)
 * n - size of the annular mask as a scalar, or vector with coordinates
 * r - radius of the annulus in pixels
 * w - width of the annulus in pixels

 #### OUTPUT
 * annulus - n x n matrix with the annulus marked with ones

 #### USAGE
```matlab
 figure;
 imshow(createAnnulus(256,32,4),[]);
```

![Basic annulus in 256x256 image with radius 32 and width 4](images/createAnnulus_basic.png)

 Create Bessel beam 2D profile
```matlab
 figure;
 imshow(log(abs(fftshift(ifft2(ifftshift(createAnnulus)))).^2+1),[]);
 colormap(gca,hot);
 caxis([0 6e-4]);
```

![Bessel pattern created by Fourier transform of default annulus created with no specified parameters](images/createAnnulus_bessel.png)

 #### REMARKS
 This could be streamlined using the bresenham circle algorithm

 Mark Kittisopikul, August 25th, 2018  
 Lab of Robert D. Goldman  
 Northwestern University

## License

    See LICENSE.txt

    Field Synthesis Demonstration - MATLAB code to demonstrate field synthesis light sheet microscopy
    Copyright (C) 2018 Reto Fioka, University of Texas Southwestern Medical Center
    Copyright (C) 2018 Mark Kittisopikul, Northwestern University

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.