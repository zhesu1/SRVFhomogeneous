# SRVFhomogeneous

Tools for geometric shape analysis of with values in homegeneous spaces based on the work by Martin Bauer, Eric Klassen and Zhe Su. 

## What is it?

This code provides tools for geometric shape analysis on open curves with values in some specific homegeneous spaces, including the 2 dimensional sphere, the 2 dimensional hyperbolic space and the space of positive definite symmetric matrices with determinant 1 (PDSM). It is able to factor out reparametrizations and rigid motions. 

For details we refer to our papers

```css
@article{SuKlBa2017,
	author = {Zhe Su, Eric Klassen, Martin Bauer},
	title = {The Square Root Velocity Framework for Curves in a Homogeneous Space},
	journal = {2017 IEEE Conference on Computer Vision and Pattern Recognition Workshops},
	pages = {680--689},
	year = 2017,
}

@article{SuKlBa2018,
  author = {Zhe Su, Eric Klassen, Martin Bauer},
	title = {Comparing curves in homogeneous spaces},
	journal = {Differential Geometry and its Applications},
	volume = 60,
	pages = {9--32},
	year = 2018
}
```

If you use our code in your work please cite our papers.

## Packages

In our code we use the dynamic programming algorithm "DynamicProgrammingQ.c" to find the optimal reparametrizations, which can be downloaded on the shape group website at FSU: [http://ssamg.stat.fsu.edu/software](http://ssamg.stat.fsu.edu/software).

The code was tested on Matlab R2018a.

## Usage

See the file "RUN_ME_sample_curves.m" in each folder of how to use the code. Compile the dynamic programming file "DynamicProgrammingQ.c" in Matlab before using programs. 

## License

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program. If not, see [http://www.gnu.org/licenses/](http://www.gnu.org/licenses/)

## Contacts

* Martin Bauer (bauer at math dot fsu dot edu)
* Zhe Su (zsu at math dot fsu dot edu)
