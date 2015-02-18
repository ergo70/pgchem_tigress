#ifndef __LIBMERCURYPP_H__
#define __LIBMERCURYPP_H__

/*
 * $Id$
 *
 * libmercury++
 *
 * A C++ library for the calculation of accurate masses
 * and abundances of isotopic peaks
 *
 * Copyright (c) 2006
 * 	Marc Kirchner <marc.kirchner@iwr.uni-heidelberg.de>
 *
 * Based on the emass implementation of Perttu Haimi
 * (see Copyright notice below).
 *
 * This code may be distributed under the terms of the
 * Lesser GNU Public License (LGPL) version 2 or any later version.
 */

/*
 *
 * Based on an algorithm developed by Alan L. Rockwood.
 *
 * Published in
 * Rockwood, A.L. and Haimi, P.: "Efficent calculation of
 * Accurate Masses of Isotopic Peaks",
 * Journal of The American Society for Mass Spectrometry
 * JASMS 03-2263, 2006
 *
 * Copyright (c) 2005 Perttu Haimi and Alan L. Rockwood
 *
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms,
 * with or without modification, are permitted provided
 * that the following conditions are met:
 *
 *    * Redistributions of source code must retain the
 *      above copyright notice, this list of conditions
 *      and the following disclaimer.
 *    * Redistributions in binary form must reproduce
 *      the above copyright notice, this list of conditions
 *      and the following disclaimer in the documentation
 *      and/or other materials provided with the distribution.
 *    * Neither the author nor the names of any contributors
 *      may be used to endorse or promote products derived
 *      from this software without specific prior written
 *      permission.
 */

/*
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND
CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,
INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR
CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT
NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#include <cstdlib>
#include <vector>
#include <cmath>

namespace mercury
{
const unsigned int MAX_ELEMENTS = 12;
const unsigned int MAX_ISOTOPES = 5;
const double electronMass = 0.00054858;
const unsigned int nIsotopes[MAX_ELEMENTS] = { 2, 2, 2, 3, 5, 2, 2, 1, 1, 1, 3 };
const double elemMasses[MAX_ELEMENTS][MAX_ISOTOPES] =
{
    {1.0078246,	2.0141021,	0,		0,	0}, // H
    {12.0, 	13.0033554, 	0, 		0,	0}, // C
    {14.0030732, 	15.0001088, 	0, 		0,	0}, // N
    {15.9949141, 	16.9991322, 	17.9991616, 	0,	0}, // O
    {31.972070, 	32.971456, 	33.967866, 	34.0,	35.967080}, // S
    {78.918336,     80.916290,  0,  0,  0}, // Br
    {34.968853,     36.965903,  0,  0,  0}, // Cl
    {126.904477,    0,  0,  0,  0}, // I
    {30.973763,     0,  0,  0,  0}, // P
    {18.998403,     0,  0,  0,  0}, // F
    {27.976928,     28.976496,  29.973772,  0,  0}, // Si
    {53.939612, 55.934939,  56.935396,  57.933278,  0} // Fe
};
const double elemAbundances[MAX_ELEMENTS][MAX_ISOTOPES] =
{
    {0.99985,	0.00015,	0,		0,	0}, // H
    {0.988930,	0.011070, 	0, 		0,	0}, // C
    {0.996337,	0.003663, 	0, 		0,	0}, // N
    {0.997590,	0.000374,	0.002036, 	0,	0}, // O
    {0.9502,	0.0075,		0.0421, 	0,	0.0002}, // S
    {0.5069,    0.4931,     0,          0,  0}, // Br
    {0.7577,    0.2423,     0,          0,  0}, // Cl
    {1.0,   0,  0,  0,  0}, // I
    {1.0,   0,  0,  0,  0}, // P
    {1.0,   0,  0,  0,  0}, // F
    {0.9223,    0.0467, 0.031,  0,  0}, // Si
    {0.058, 0.9172, 0.022,  0.0028, 0} // Fe
};

/*
	mercury:	calculates the expected isotpic distribution
			for a given composition
	parameters:
	msa_mz		returned mass list
	msa_abundance	returned abundance list
	composition	a vector of length MAX_ELEMENTS giving the
			the number for occurences of each element
	charge		the charge state for which the isotopic pattern
			is to be calculated
	limit		a pruning limit, (Rockwood et al. 2006) suggest
			a value between 10e-25 and 10e-30
*/
int mercury(std::vector<double>& msa_mz, std::vector<double>& msa_abundance, const std::vector<unsigned int>& composition, const int charge, const double limit);
}
#endif
