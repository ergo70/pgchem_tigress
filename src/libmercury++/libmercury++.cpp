#include "libmercury++.h"

// declare private stuff
namespace mercury {
	void convolve(std::vector<double>& result_mz, std::vector<double>& result_ab, const std::vector<double>& mz1, const std::vector<double>& ab1, const std::vector<double>& mz2, const std::vector<double>& ab2);
	void prune(std::vector<double>& mz, std::vector<double>& ab, const double limit);
}


void mercury::convolve(std::vector<double>& result_mz, std::vector<double>& result_ab, const std::vector<double>& mz1, const std::vector<double>& ab1, const std::vector<double>& mz2, const std::vector<double>& ab2)
{
	size_t n1 = mz1.size();
	size_t n2 = mz2.size();
	if ((n1+n2) == 0)
		return;

	result_mz.clear();
	result_ab.clear();
	// the following two lines speed up calculations by a factor of 3 (!)
	result_mz.resize(n1+n2);
	result_ab.resize(n1+n2);
	// for each isotope peak in the compound...
	double totalAbundance, massExpectation, ithMass, ithAbundance;
	for (size_t k = 0; k < n1 + n2 - 1; k++) {
		totalAbundance = 0;
		massExpectation = 0;
		size_t start = k < (n2 - 1) ? 0 : k - n2 + 1; // start=max(0, k-n2+1)
		size_t end = k < (n1 - 1) ? k : n1 - 1;       // end=min(n1 - 1, k)
		// ... calculate the mass expectation value and the abundance
		for (size_t i = start; i <= end; i++) {
			ithAbundance = ab1[i] * ab2[k - i];
			if (ithAbundance > 0) {
				totalAbundance += ithAbundance;
				ithMass = mz1[i] + mz2[k - i];
				massExpectation += ithAbundance * ithMass;
			}
		}
		// do NOT throw away isotopes with zero probability, this would
		// screw up the isotope count k !!
		result_mz[k] = totalAbundance > 0 ? massExpectation / totalAbundance : 0;
		result_ab[k] = totalAbundance;
	}
}

void mercury::prune(std::vector<double>& mz, std::vector<double>& ab, const double limit)
{
	size_t i;
	for (i = 0; i < ab.size(); i++) {
		if(ab[i] > limit)
			break;
	}
	mz.erase(mz.begin(), mz.begin()+i);
	ab.erase(ab.begin(), ab.begin()+i);

	// prune the end
	for (i = ab.size()-1; i >= 0; i--) {
		if(ab[i] > limit)
			break;
	}
	mz.resize(i+1);
	ab.resize(i+1);
	//mz.erase(mz.begin()+i+1, mz.end());
	//ab.erase(ab.begin()+i+1, ab.end());
}

int mercury::mercury(std::vector<double>& msa_mz, std::vector<double>& msa_abundance, const std::vector<unsigned int>& composition, const int charge, const double limit)
{
	if (composition.size() != MAX_ELEMENTS) {
		return(-1);
	}

	unsigned int n;
	std::vector<double> tmp_mz, tmp_abundance, esa_mz, esa_abundance ;
	bool msa_initialized = false;

	// walk through the elements
	for (unsigned int e = 0; e < MAX_ELEMENTS; e++) {
		// if the element is present in the composition,
		// then calculate ESA and update MSA
//		std::cout << "e=" << e << std::endl;
		n = composition[e];
		if (n) {
			// initialize ESA
			esa_mz.assign(elemMasses[e], elemMasses[e]+nIsotopes[e]);
			esa_abundance.assign(elemAbundances[e], elemAbundances[e]+nIsotopes[e]);
		//	esa_mz.resize(n*nIsotopes[e]);
		//	esa_abundance.resize(n*nIsotopes[e]);
			//while (n > 0) {
			while (1) {
//				std::cout << "n=" << n << std::endl;
				// check if we need to do the MSA update
				if (n & 1) {
					// MSA update
					if (msa_initialized) {
						// normal update
						convolve(tmp_mz, tmp_abundance, msa_mz, msa_abundance, esa_mz, esa_abundance);
						msa_mz = tmp_mz;
						msa_abundance = tmp_abundance;
					} else {
						// for the first assignment MSA=ESA
						msa_mz = esa_mz;
						msa_abundance = esa_abundance;
						msa_initialized = true;
					}
					prune(msa_mz, msa_abundance, limit);
				}
				// the ESA update is always carried out (with the exception of
				// the last time, i.e. when n==1)
				if (n==1) {
					break;
				}
				convolve(tmp_mz, tmp_abundance, esa_mz, esa_abundance, esa_mz, esa_abundance);
				esa_mz = tmp_mz;
				esa_abundance = tmp_abundance;
				prune(esa_mz, esa_abundance, limit);
				n = n >> 1;
			}
		}
	}
	// take charge into account (placing the if around two loops is faster
	// than vice versa
	if (charge > 0) {
		for (std::vector<double>::iterator i = msa_mz.begin(); i != msa_mz.end(); i++) {
			*i = *i / abs(charge) - electronMass;
		}
	}
	if (charge < 0) {
		for (std::vector<double>::iterator i = msa_mz.begin(); i != msa_mz.end(); i++) {
			*i = *i / abs(charge) + electronMass;
		}
	}
	return(0);
}
