//////////////////////////
// tools.hpp
//////////////////////////
//
// Collection of analysis tools for gevolution
//
// Author: Julian Adamek (Université de Genève & Observatoire de Paris & Queen Mary University of London)
//
// Last modified: April 2019
//
//////////////////////////

#ifndef TOOLS_HEADER
#define TOOLS_HEADER

#ifndef Cplx
#define Cplx Imag
#endif

#define KTYPE_GRID      0
#define KTYPE_LINEAR    1

using namespace std;
using namespace LATfield2;


#ifdef FFT3D
//////////////////////////
// extractCrossSpectrum
//////////////////////////
// Description:
//   generates the cross spectrum for two Fourier images
//
// Arguments:
//   fld1FT     reference to the first Fourier image for which the cross spectrum should be extracted
//   fld2FT     reference to the second Fourier image for which the cross spectrum should be extracted
//   kbin       allocated array that will contain the central k-value for the bins
//   power      allocated array that will contain the average power in each bin
//   kscatter   allocated array that will contain the k-scatter for each bin
//   pscatter   allocated array that will contain the scatter in power for each bin
//   occupation allocated array that will count the number of grid points contributing to each bin
//   numbin     number of bins (minimum size of all arrays)
//   ktype      flag indicating which definition of momentum to be used
//                  0: grid momentum
//                  1: linear (default)
//   comp1      for component-wise cross spectra, the component for the first field (ignored if negative)
//   comp2      for component-wise cross spectra, the component for the second field (ignored if negative)
//
// Returns:
//
//////////////////////////
template<int mu_bins = 1>
void extractCrossSpectrum(Field<Cplx> & fld1FT, Field<Cplx> & fld2FT, Real * kbin, Real * power, Real * kscatter, Real * pscatter, int * occupation, const int numbins, const bool deconvolve = true, const int ktype = KTYPE_LINEAR, const int comp1 = -1, const int comp2 = -1)
{
	int i, weight, j;//j index of mu_bins
	const int linesize = fld1FT.lattice().size(1);
	Real * typek2;
	Real * sinc;
	Real k2max, k2, s;
	rKSite k(fld1FT.lattice());
	Cplx p;

	typek2 = (Real *) malloc(linesize * sizeof(Real));
	sinc = (Real *) malloc(linesize * sizeof(Real));

	if (ktype == KTYPE_GRID)
	{
		for (i = 0; i < linesize; i++)
		{
			typek2[i] = 2. * (Real) linesize * sin(M_PI * (Real) i / (Real) linesize);
			typek2[i] *= typek2[i];
		}
	}
	else
	{
		for (i = 0; i <= linesize/2; i++)
		{
			typek2[i] = 2. * M_PI * (Real) i;
			typek2[i] *= typek2[i];
		}
		for (; i < linesize; i++)
		{
			typek2[i] = 2. * M_PI * (Real) (linesize-i);
			typek2[i] *= typek2[i];
		}
	}

	sinc[0] = 1.;
	if (deconvolve)
	{
		for (i = 1; i <= linesize / 2; i++)
		{
			sinc[i] = sin(M_PI * (float) i / (float) linesize) * (float) linesize / (M_PI * (float) i);
		}
	}
	else
	{
		for (i = 1; i <= linesize / 2; i++)
		{
			sinc[i] = 1.;
		}
	}
	for (; i < linesize; i++)
	{
		sinc[i] = sinc[linesize-i];
	}

	k2max = 3. * typek2[linesize/2];

	for (i = 0; i < numbins*mu_bins; i++)
	{
		kbin[i] = 0.;
		power[i] = 0.;
		kscatter[i] = 0.;
		pscatter[i] = 0.;
		occupation[i] = 0;
	}

	for (k.first(); k.test(); k.next())
	{
		if (k.coord(0) == 0 && k.coord(1) == 0 && k.coord(2) == 0)
			continue;
		else if (k.coord(0) == 0)
			weight = 1;
		else if ((k.coord(0) == linesize/2) && (linesize % 2 == 0))
			weight = 1;
		else
			weight = 2;

		k2 = typek2[k.coord(0)] + typek2[k.coord(1)] + typek2[k.coord(2)];
		s = sinc[k.coord(0)] * sinc[k.coord(1)] * sinc[k.coord(2)]; // What is s and sinc ?
		s *= s;

		if (comp1 >= 0 && comp2 >= 0 && comp1 < fld1FT.components() && comp2 < fld2FT.components())
		{
			p = fld1FT(k, comp1) * fld2FT(k, comp2).conj();
		}
		else if (fld1FT.symmetry() == LATfield2::symmetric)
		{
			p = fld1FT(k, 0, 1) * fld2FT(k, 0, 1).conj();
			p += fld1FT(k, 0, 2) * fld2FT(k, 0, 2).conj();
			p += fld1FT(k, 1, 2) * fld2FT(k, 1, 2).conj();
			p *= 2.;
			p += fld1FT(k, 0, 0) * fld2FT(k, 0, 0).conj();
			p += fld1FT(k, 1, 1) * fld2FT(k, 1, 1).conj();
			p += fld1FT(k, 2, 2) * fld2FT(k, 2, 2).conj();
		}
		else
		{
			p = Cplx(0., 0.);
			for (i = 0; i < fld1FT.components(); i++)
				p += fld1FT(k, i) * fld2FT(k, i).conj();
		}

		i = (int) floor((double) ((Real) numbins * sqrt(k2 / k2max)));
		if (i < numbins) //ther is no else except k = kmax
		{
			j = (int) floor((double) ((Real) (mu_bins) * sqrt( typek2[k.coord(0)] / k2)));
			if(j >= mu_bins) j = mu_bins-1; //It can't be bigger?
			kbin[i * mu_bins + j] += weight * sqrt(k2);
			kscatter[i * mu_bins + j] += weight * k2;
			power[i * mu_bins + j] += weight * p.real() * k2 * sqrt(k2) / s;//ask Julien what is s? and for multipole, should we remove k2^3/2 from here? or is there other way to do it later?
			pscatter[i * mu_bins + j] += weight * p.real() * p.real() * k2 * k2 * k2 / s / s;
			occupation[i * mu_bins + j] += weight;
		}
	}

	free(typek2);
	free(sinc);

	if (parallel.isRoot())
	{
		// Here different cpus communicate, each of them does the same calculation for different part of the lattice
#ifdef SINGLE
		MPI_Reduce(MPI_IN_PLACE, (void *) kbin, numbins*mu_bins, MPI_FLOAT, MPI_SUM, 0, parallel.lat_world_comm());
		MPI_Reduce(MPI_IN_PLACE, (void *) kscatter, numbins*mu_bins, MPI_FLOAT, MPI_SUM, 0, parallel.lat_world_comm());
		MPI_Reduce(MPI_IN_PLACE, (void *) power, numbins*mu_bins, MPI_FLOAT, MPI_SUM, 0, parallel.lat_world_comm());
		MPI_Reduce(MPI_IN_PLACE, (void *) pscatter, numbins*mu_bins, MPI_FLOAT, MPI_SUM, 0, parallel.lat_world_comm());
#else
		MPI_Reduce(MPI_IN_PLACE, (void *) kbin, numbins*mu_bins, MPI_DOUBLE, MPI_SUM, 0, parallel.lat_world_comm());
		MPI_Reduce(MPI_IN_PLACE, (void *) kscatter, numbins*mu_bins, MPI_DOUBLE, MPI_SUM, 0, parallel.lat_world_comm());
		MPI_Reduce(MPI_IN_PLACE, (void *) power, numbins*mu_bins, MPI_DOUBLE, MPI_SUM, 0, parallel.lat_world_comm());
		MPI_Reduce(MPI_IN_PLACE, (void *) pscatter, numbins*mu_bins, MPI_DOUBLE, MPI_SUM, 0, parallel.lat_world_comm());
#endif
		MPI_Reduce(MPI_IN_PLACE, (void *) occupation, numbins*mu_bins, MPI_INT, MPI_SUM, 0, parallel.lat_world_comm());
		//Normal p(k)
		for (i = 0; i < numbins*mu_bins; i++)
		{
			if (occupation[i] > 0)
			{
				kscatter[i] = sqrt(kscatter[i] * occupation[i] - kbin[i] * kbin[i]) / occupation[i];
				if (!isfinite(kscatter[i])) kscatter[i] = 0.;
				kbin[i] = kbin[i] / occupation[i];
				power[i] /= occupation[i];
				pscatter[i] = sqrt(pscatter[i] / occupation[i] - power[i] * power[i]);
				if (!isfinite(pscatter[i])) pscatter[i] = 0.;
			}
		}
		//here do the for loop for multipole calculations, integrals for multipole analytically(trapz, begining and the end point should be done rectangular), we should go forward in mu-bin until it is not zero
		//-------
		if(mu_bins>1){
		for (i = 0; i < numbins; i++) // this fixes k
		{
			Real mu_minus, mu_plus, P_hexadecapole, P_monopole, P_quadrupole;
			int m=0;//to count how many bins to skip, m<mu_bins
			while(m<mu_bins && occupation[i*mu_bins+m] == 0 ) m++;
			P_monopole = power[i*mu_bins + m] *((m+0.5)/mu_bins);//assuming constant value up to first non-empty
			//((1 + 2 m) (1 + 4 m + 4 m^2 - 4 mubins^2))/(16 mubins^3)
			P_quadrupole = power[i*mu_bins + m ]*(1+2*m)*(1 + 4*m + 4*m*m - 4*mu_bins*mu_bins )/(16*pow(mu_bins,3) );
			//(14 (1/2 + m)^5 - 20 (1/2 + m)^3 mubins^2 + 3 mubins^4 + 6 m mubins^4)/(16 mubins^5)
			P_hexadecapole = power[i*mu_bins + m ]* (14*pow((0.5 + m),5) - 20* pow((0.5 + m),3)*mu_bins*mu_bins + 3*pow(mu_bins,4) + 6*m*pow(mu_bins,4) )
			/(16*pow(mu_bins,5));
			// for quad and hexa we should integrate Legendre (from Mathematica)
			//do the same for hexa
			//for (j=0, j<5,j++)
			for(j=m; true ;j=m)
			{
				m ++;
				while(m<mu_bins && occupation[i*mu_bins+m]==0 ) m++; //next empty
				if(m >= mu_bins) break;//
				mu_minus = (1./mu_bins)*(j+0.5);//previous non-empty bin
				mu_plus = (1./mu_bins)*(m +0.5);//current non-empty bin
				P_monopole += (power[ i * mu_bins + j] + power[i * mu_bins + m])*(mu_plus - mu_minus)/2.;
				//write powers either using pow() function or use *
				P_quadrupole += (mu_plus - mu_minus)*( (-2 + 3*(mu_minus)*(mu_minus) + 2*mu_minus*mu_plus + mu_plus*mu_plus )*power[ i * mu_bins + j]
				+ (-2 + mu_minus*mu_minus + 2*mu_minus*mu_plus + 3*mu_plus*mu_plus)*power[i * mu_bins + m] )/8.;
				//
				P_hexadecapole += (mu_plus - mu_minus)*( (9 + 35*pow(mu_minus,4) + 28*pow(mu_minus,3)*mu_plus - 15*mu_plus*mu_plus
					 + 7*pow(mu_plus,4) +3*mu_minus*mu_minus*(-15 + 7*mu_plus*mu_plus) +2*mu_minus*mu_plus*(-15 + 7*mu_plus*mu_plus))*power[ i * mu_bins + j]
					 + (9 + 7*pow(mu_minus,4) +14*pow(mu_minus,3)*mu_plus - 45*mu_plus*mu_plus + 35*pow(mu_plus,4) +3*mu_minus*mu_minus*(-5 + 7*mu_plus*mu_plus)
					 +mu_minus*(-30*mu_plus + 28*pow(mu_plus,3) ))*power[i * mu_bins + m] )/48.;
			}
			//adding the the last bin, ask Julian, why should it be j which is the first non-zero bin?
			P_monopole += (power[i*mu_bins + j] )*(1-(j+0.5)/mu_bins);
			//quad: -(((1 + 2 m) (1 + 4 m + 4 m^2 - 4 mubins^2))/(16 mubins^3))
			P_quadrupole += -power[i*mu_bins + j]*(1+2*j)*(1+4*j + 4*j*j - 4*mu_bins*mu_bins)/(16*pow(mu_bins,3));
			//P_hexadecapole: (-7 (1 + 2 m)^5 + 40 (1 + 2 m)^3 mubins^2 - 48 (1 + 2 m) mubins^4)/(256 mubins^5)
			P_hexadecapole += power[i*mu_bins + j]*(-7*pow((1 + 2*j),5) + 40*pow((1 + 2*j),3)*mu_bins*mu_bins - 48*(1 + 2*j)*pow(mu_bins,4))/(256*pow(mu_bins,5) );
			//do the same for quad and hexa
			for(j=1;j<mu_bins;j++) {
				//
				if(occupation [i*mu_bins + j] > 0) kbin[i*mu_bins] = (kbin[i*mu_bins]*occupation[i*mu_bins] + kbin[i*mu_bins+j]*occupation[i*mu_bins+j])/(occupation[i*mu_bins]+occupation[i*mu_bins+j]);
				if(occupation [i*mu_bins + j] > 0) kscatter[i*mu_bins] = (kscatter[i*mu_bins]*occupation[i*mu_bins] + kscatter[i*mu_bins+j]*occupation[i*mu_bins+j])/(occupation[i*mu_bins]+occupation[i*mu_bins+j]);

				occupation[i*mu_bins]+= occupation[i*mu_bins + j];
			}
			//we are
			power[i*mu_bins+0] = P_monopole;
			power[i*mu_bins+1] = P_quadrupole;
			power[i*mu_bins+2] = P_hexadecapole;
		}
	}
		//-------
	}
	else
	{
#ifdef SINGLE
		MPI_Reduce((void *) kbin, NULL, numbins*mu_bins, MPI_FLOAT, MPI_SUM, 0, parallel.lat_world_comm());
		MPI_Reduce((void *) kscatter, NULL, numbins*mu_bins, MPI_FLOAT, MPI_SUM, 0, parallel.lat_world_comm());
		MPI_Reduce((void *) power, NULL, numbins*mu_bins, MPI_FLOAT, MPI_SUM, 0, parallel.lat_world_comm());
		MPI_Reduce((void *) pscatter, NULL, numbins*mu_bins, MPI_FLOAT, MPI_SUM, 0, parallel.lat_world_comm());
#else
		MPI_Reduce((void *) kbin, NULL, numbins*mu_bins, MPI_DOUBLE, MPI_SUM, 0, parallel.lat_world_comm());
		MPI_Reduce((void *) kscatter, NULL, numbins*mu_bins, MPI_DOUBLE, MPI_SUM, 0, parallel.lat_world_comm());
		MPI_Reduce((void *) power, NULL, numbins*mu_bins, MPI_DOUBLE, MPI_SUM, 0, parallel.lat_world_comm());
		MPI_Reduce((void *) pscatter, NULL, numbins*mu_bins, MPI_DOUBLE, MPI_SUM, 0, parallel.lat_world_comm());
#endif
		MPI_Reduce((void *) occupation, NULL, numbins*mu_bins, MPI_INT, MPI_SUM, 0, parallel.lat_world_comm());
	}
}


//////////////////////////
// extractPowerSpectrum
//////////////////////////
// Description:
//   generates the power spectrum for a Fourier image
//
// Arguments:
//   fldFT      reference to the Fourier image for which the power spectrum should be extracted
//   kbin       allocated array that will contain the central k-value for the bins
//   power      allocated array that will contain the average power in each bin
//   kscatter   allocated array that will contain the k-scatter for each bin
//   pscatter   allocated array that will contain the scatter in power for each bin
//   occupation allocated array that will count the number of grid points contributing to each bin
//   numbin     number of bins (minimum size of all arrays)
//   ktype      flag indicating which definition of momentum to be used
//                  0: grid momentum
//                  1: linear (default)
//
// Returns:
//
//////////////////////////
template<int mu_bins = 1>
void extractPowerSpectrum(Field<Cplx> & fldFT, Real * kbin, Real * power, Real * kscatter, Real * pscatter, int * occupation, const int numbins, const bool deconvolve = true, const int ktype = KTYPE_LINEAR)
{
	extractCrossSpectrum<mu_bins>(fldFT, fldFT, kbin, power, kscatter, pscatter, occupation, numbins, deconvolve, ktype);
}
#endif


//////////////////////////
// writePowerSpectrum
//////////////////////////
// Description:
//   writes power spectra as tabulated data into ASCII file
//
// Arguments:
//   kbin           array containing the central values of k for each bin
//   power          array containing the central values of P(k) for each bin
//   kscatter       array containing the statistical error on k for each bin
//   pscatter       array containing the statistical error on P(k) for each bin
//   occupation     array containing the number of k-modes contributing to each bin
//   numbins        total number of bins (length of the arrays)
//   rescalek       unit conversion factor for k
//   rescalep       unit conversion factor for P(k)
//   filename       output file name
//   description    descriptive header
//   a              scale factor for this spectrum
//   z_target       target redshift for this output (used only if EXACT_OUTPUT_REDSHIFTS is defined)
//
// Returns:
//
//////////////////////////
//changed
template<int mu_bins = 1>
void writePowerSpectrum(Real * kbin, Real * power, Real * kscatter, Real * pscatter, int * occupation, const int numbins, const Real rescalek, const Real rescalep, const char * filename, const char * description, double a, const double z_target = -1)
{
	if (parallel.isRoot())
	{
		//interpolate between two time step that the ... (read the previous outputted file)
#ifdef EXACT_OUTPUT_REDSHIFTS
//changed
		Real * power2 = (Real *) malloc(numbins * sizeof(Real)*mu_bins );

		for (int i = 0; i < numbins*mu_bins; i++)
			power2[i] = power[i]/rescalep;

		if (1. / a < z_target + 1.)
		{
			FILE * infile = fopen(filename, "r");
			double weight = 1.;
			int count = 0;
			if (infile != NULL)
			{
				fscanf(infile, "%*[^\n]\n");
				if (fscanf(infile, "# redshift z=%lf\n", &weight) != 1)
				{
					cout << " error parsing power spectrum file header for interpolation (EXACT_OUTPUT_REDSHIFTS)" << endl;
					weight = 1.;
				}
				else
				{
					weight = (weight - z_target) / (1. + weight - 1./a);
					fscanf(infile, "%*[^\n]\n");
					//changed
					for (int i = 0; i < numbins*mu_bins; i+=mu_bins)
					{
						if (occupation[i] > 0)
						{
//changed
if(mu_bins>1){
	#ifdef SINGLE
								if(fscanf(infile, " %*e %*e %e %e %e %*d \n", power2+i, power2+i+1, power2+i+2) != 3)
	#else
								if(fscanf(infile, " %*e %*e %le %le %le %*d \n", power2+i, power2+i+1, power2+i+2) != 3)
	#endif
	{
		cout << " error parsing power spectrum file data " << i << " for interpolation (EXACT_OUTPUT_REDSHIFTS)" << endl;
		break;
	}
	else count++;
}
else{
	#ifdef SINGLE
								if(fscanf(infile, " %*e %e %*e %*e %*d \n", power2+i) != 1)
	#else
								if(fscanf(infile, " %*e %le %*e %*e %*d \n", power2+i) != 1)
	#endif
	{
		cout << " error parsing power spectrum file data " << i << " for interpolation (EXACT_OUTPUT_REDSHIFTS)" << endl;
		break;
	}
	else count++;
}
						}
					}
				}
				fclose(infile);
					//changed
				for (int i = 0; i < numbins*mu_bins; i++)
					power2[i] = (1.-weight)*power2[i] + weight*power[i]/rescalep;

				a = 1. / (z_target + 1.);
			}
		}
#endif // EXACT_OUTPUT_REDSHIFTS
		FILE * outfile = fopen(filename, "w");
		if (outfile == NULL)
		{
			cout << " error opening file for power spectrum output!" << endl;
		}
		else
		{
			fprintf(outfile, "# %s\n", description);
			fprintf(outfile, "# redshift z=%f\n", (1./a)-1.);
			if(mu_bins>1) fprintf(outfile, "# k              sigma(k)              P0              P2              P4              count\n");
			else fprintf(outfile, "# k              Pk             sigma(k)       sigma(Pk)      count\n");
			//changed
			for (int i = 0; i < numbins; i++)
			{
				if (occupation[i*mu_bins] > 0){
						#ifdef EXACT_OUTPUT_REDSHIFTS
											//changed
											if(mu_bins>1) fprintf(outfile, "  %e   %e   %e   %e   %e   %d\n", kbin[i*mu_bins]/rescalek, kscatter[i*mu_bins]/rescalek,power2[i*mu_bins],power2[i*mu_bins+1],power2[i*mu_bins+2], occupation[i*mu_bins]);
											else fprintf(outfile, "  %e   %e   %e   %e   %d\n", kbin[i]/rescalek, power2[i], kscatter[i]/rescalek, pscatter[i]/rescalep/ sqrt(occupation[i]), occupation[i]);
						#else
											if (mu_bins>1) fprintf(outfile, "  %e   %e   %e   %e   %e   %d\n", kbin[i*mu_bins]/rescalek, kscatter[i*mu_bins]/rescalek, power[i*mu_bins]/rescalep, power[i*mu_bins+1]/rescalep,power[i*mu_bins+2]/rescalep, occupation[i*mu_bins]);
											else fprintf(outfile, "  %e   %e   %e   %e   %d\n", kbin[i]/rescalek, power[i]/rescalep, kscatter[i]/rescalek, pscatter[i]/rescalep/ sqrt(occupation[i]), occupation[i]);
						#endif
						}
			}
			fclose(outfile);
		}
#ifdef EXACT_OUTPUT_REDSHIFTS
		free(power2);
#endif
	}
}


//////////////////////////
// computeVectorDiagnostics
//////////////////////////
// Description:
//   computes some diagnostics for the spin-1 perturbation
//
// Arguments:
//   Bi         reference to the real-space vector field to analyze
//   mdivB      will contain the maximum value of the divergence of Bi
//   mcurlB     will contain the maximum value of the curl of Bi
//
// Returns:
//
//////////////////////////

void computeVectorDiagnostics(Field<Real> & Bi, Real & mdivB, Real & mcurlB)
{
	Real b1, b2, b3, b4;
	const Real linesize = (Real) Bi.lattice().sizeLocal(0);
	Site x(Bi.lattice());

	mdivB = 0.;
	mcurlB = 0.;

	for (x.first(); x.test(); x.next())
	{
		b1 = fabs((Bi(x,0)-Bi(x-0,0)) + (Bi(x,1)-Bi(x-1,1)) + (Bi(x,2)-Bi(x-2,2))) * linesize;
		if (b1 > mdivB) mdivB = b1;
		b1 = 0.5 * (Bi(x,0) + Bi(x+0,1) - Bi(x+1,0) - Bi(x,1) + Bi(x+2,0) + Bi(x+0+2,1) - Bi(x+1+2,0) - Bi(x+2,1)) * linesize;
		b2 = 0.5 * (Bi(x,0) + Bi(x+0,2) - Bi(x+2,0) - Bi(x,2) + Bi(x+1,0) + Bi(x+0+1,2) - Bi(x+2+1,0) - Bi(x+1,2)) * linesize;
		b3 = 0.5 * (Bi(x,2) + Bi(x+2,1) - Bi(x+1,2) - Bi(x,1) + Bi(x+0,2) + Bi(x+2+0,1) - Bi(x+1+0,2) - Bi(x+0,1)) * linesize;
		b4 = sqrt(b1 * b1 + b2 * b2 + b3 * b3);
		if (b4 > mcurlB) mcurlB = b4;
	}

	parallel.max<Real>(mdivB);
	parallel.max<Real>(mcurlB);
}


//////////////////////////
// computeTensorDiagnostics
//////////////////////////
// Description:
//   computes some diagnostics for the spin-2 perturbation
//
// Arguments:
//   hij        reference to the real-space tensor field to analyze
//   mdivh      will contain the maximum value of the divergence of hij
//   mtraceh    will contain the maximum value of the trace of hij
//   mnormh     will contain the maximum value of the norm of hij
//
// Returns:
//
//////////////////////////

void computeTensorDiagnostics(Field<Real> & hij, Real & mdivh, Real & mtraceh, Real & mnormh)
{
	Real d1, d2, d3;
	const Real linesize = (Real) hij.lattice().sizeLocal(0);
	Site x(hij.lattice());

	mdivh = 0.;
	mtraceh = 0.;
	mnormh = 0.;

	for (x.first(); x.test(); x.next())
	{
		d1 = (hij(x+0, 0, 0) - hij(x, 0, 0) + hij(x, 0, 1) - hij(x-1, 0, 1) + hij(x, 0, 2) - hij(x-2, 0, 2)) * linesize;
		d2 = (hij(x+1, 1, 1) - hij(x, 1, 1) + hij(x, 0, 1) - hij(x-0, 0, 1) + hij(x, 1, 2) - hij(x-2, 1, 2)) * linesize;
		d3 = (hij(x+2, 2, 2) - hij(x, 2, 2) + hij(x, 0, 2) - hij(x-0, 0, 2) + hij(x, 1, 2) - hij(x-1, 1, 2)) * linesize;
		d1 = sqrt(d1 * d1 + d2 * d2 + d3 * d3);
		if (d1 > mdivh) mdivh = d1;
		d1 = fabs(hij(x, 0, 0) + hij(x, 1, 1) + hij(x, 2, 2));
		if (d1 > mtraceh) mtraceh = d1;
		d1 = sqrt(hij(x, 0, 0) * hij(x, 0, 0) + 2. * hij(x, 0, 1) * hij(x, 0, 1) + 2. * hij(x, 0, 2)* hij(x, 0, 2) + hij(x, 1, 1) * hij(x, 1, 1) + 2. * hij(x, 1, 2) * hij(x, 1, 2) + hij(x, 2, 2) * hij(x, 2, 2));
		if (d1 > mnormh) mnormh = d1;
	}

	parallel.max<Real>(mdivh);
	parallel.max<Real>(mtraceh);
	parallel.max<Real>(mnormh);
}


//////////////////////////
// findIntersectingLightcones
//////////////////////////
// Description:
//   determines periodic copies of light cone vertex for which the present
//   look-back interval may overlap with a given spatial domain
//
// Arguments:
//   lightcone  reference to structure describing light cone geometry
//   outer      outer (far) limit of look-back interval
//   inner      inner (close) limit of look-back interval
//   domain     array of domain boundaries
//   vertex     will contain array of relevant vertex locations
//
// Returns:
//   number of vertices found
//
//////////////////////////

int findIntersectingLightcones(lightcone_geometry & lightcone, double outer, double inner, double * domain, double vertex[MAX_INTERSECTS][3])
{
	int range = (int) ceil(outer) + 1;
	int u, v, w, n = 0;
	double corner[8][3];
	double rdom, dist;

	corner[0][0] = domain[0];
	corner[0][1] = domain[1];
	corner[0][2] = domain[2];

	corner[1][0] = domain[3];
	corner[1][1] = domain[1];
	corner[1][2] = domain[2];

	corner[2][0] = domain[0];
	corner[2][1] = domain[4];
	corner[2][2] = domain[2];

	corner[3][0] = domain[3];
	corner[3][1] = domain[4];
	corner[3][2] = domain[2];

	corner[4][0] = domain[0];
	corner[4][1] = domain[1];
	corner[4][2] = domain[5];

	corner[5][0] = domain[3];
	corner[5][1] = domain[1];
	corner[5][2] = domain[5];

	corner[6][0] = domain[0];
	corner[6][1] = domain[4];
	corner[6][2] = domain[5];

	corner[7][0] = domain[3];
	corner[7][1] = domain[4];
	corner[7][2] = domain[5];

	for (u = -range; u <= range; u++)
	{
		for (v = -range; v <= range; v++)
		{
			for (w = -range; w <= range; w++)
			{
				if (n >= MAX_INTERSECTS)
				{
					cout << COLORTEXT_YELLOW << " /!\\ warning" << COLORTEXT_RESET << ": maximum number of lightcone intersects exceeds MAX_INTERSECTS = " << MAX_INTERSECTS << " for domain (" << domain[0] << ", " << domain[1] << ", " << domain[2] << ") - (" << domain[3] << ", " << domain[4] << ", " << domain[5] << "); some data may be missing in output!" << endl;
					return MAX_INTERSECTS;
				}
				vertex[n][0] = lightcone.vertex[0] + u;
				vertex[n][1] = lightcone.vertex[1] + v;
				vertex[n][2] = lightcone.vertex[2] + w;

				// first, check if domain lies outside outer sphere
				if (vertex[n][0] < domain[0])
				{
					if (vertex[n][1] < domain[1])
					{
						if (vertex[n][2] < domain[2])
						{
							if (sqrt((vertex[n][0]-corner[0][0])*(vertex[n][0]-corner[0][0]) + (vertex[n][1]-corner[0][1])*(vertex[n][1]-corner[0][1]) + (vertex[n][2]-corner[0][2])*(vertex[n][2]-corner[0][2])) > outer) continue;
						}
						else if (vertex[n][2] > domain[5])
						{
							if (sqrt((vertex[n][0]-corner[4][0])*(vertex[n][0]-corner[4][0]) + (vertex[n][1]-corner[4][1])*(vertex[n][1]-corner[4][1]) + (vertex[n][2]-corner[4][2])*(vertex[n][2]-corner[4][2])) > outer) continue;
						}
						else if (sqrt((vertex[n][0]-domain[0])*(vertex[n][0]-domain[0]) + (vertex[n][1]-domain[1])*(vertex[n][1]-domain[1])) > outer) continue;
					}
					else if (vertex[n][1] > domain[4])
					{
						if (vertex[n][2] < domain[2])
						{
							if (sqrt((vertex[n][0]-corner[2][0])*(vertex[n][0]-corner[2][0]) + (vertex[n][1]-corner[2][1])*(vertex[n][1]-corner[2][1]) + (vertex[n][2]-corner[2][2])*(vertex[n][2]-corner[2][2])) > outer) continue;
						}
						else if (vertex[n][2] > domain[5])
						{
							if (sqrt((vertex[n][0]-corner[6][0])*(vertex[n][0]-corner[6][0]) + (vertex[n][1]-corner[6][1])*(vertex[n][1]-corner[6][1]) + (vertex[n][2]-corner[6][2])*(vertex[n][2]-corner[6][2])) > outer) continue;
						}
						else if (sqrt((vertex[n][0]-domain[0])*(vertex[n][0]-domain[0]) + (vertex[n][1]-domain[4])*(vertex[n][1]-domain[4])) > outer) continue;
					}
					else
					{
						if (vertex[n][2] < domain[2])
						{
							if (sqrt((vertex[n][0]-domain[0])*(vertex[n][0]-domain[0]) + (vertex[n][2]-domain[2])*(vertex[n][2]-domain[2])) > outer) continue;
						}
						else if (vertex[n][2] > domain[5])
						{
							if (sqrt((vertex[n][0]-domain[0])*(vertex[n][0]-domain[0]) + (vertex[n][2]-domain[5])*(vertex[n][2]-domain[5])) > outer) continue;
						}
						else if (domain[0]-vertex[n][0] > outer) continue;
					}
				}
				else if (vertex[n][0] > domain[3])
				{
					if (vertex[n][1] < domain[1])
					{
						if (vertex[n][2] < domain[2])
						{
							if (sqrt((vertex[n][0]-corner[1][0])*(vertex[n][0]-corner[1][0]) + (vertex[n][1]-corner[1][1])*(vertex[n][1]-corner[1][1]) + (vertex[n][2]-corner[1][2])*(vertex[n][2]-corner[1][2])) > outer) continue;
						}
						else if (vertex[n][2] > domain[5])
						{
							if (sqrt((vertex[n][0]-corner[5][0])*(vertex[n][0]-corner[5][0]) + (vertex[n][1]-corner[5][1])*(vertex[n][1]-corner[5][1]) + (vertex[n][2]-corner[5][2])*(vertex[n][2]-corner[5][2])) > outer) continue;
						}
						else if (sqrt((vertex[n][0]-domain[3])*(vertex[n][0]-domain[3]) + (vertex[n][1]-domain[1])*(vertex[n][1]-domain[1])) > outer) continue;
					}
					else if (vertex[n][1] > domain[4])
					{
						if (vertex[n][2] < domain[2])
						{
							if (sqrt((vertex[n][0]-corner[3][0])*(vertex[n][0]-corner[3][0]) + (vertex[n][1]-corner[3][1])*(vertex[n][1]-corner[3][1]) + (vertex[n][2]-corner[3][2])*(vertex[n][2]-corner[3][2])) > outer) continue;
						}
						else if (vertex[n][2] > domain[5])
						{
							if (sqrt((vertex[n][0]-corner[7][0])*(vertex[n][0]-corner[7][0]) + (vertex[n][1]-corner[7][1])*(vertex[n][1]-corner[7][1]) + (vertex[n][2]-corner[7][2])*(vertex[n][2]-corner[7][2])) > outer) continue;
						}
						else if (sqrt((vertex[n][0]-domain[3])*(vertex[n][0]-domain[3]) + (vertex[n][1]-domain[4])*(vertex[n][1]-domain[4])) > outer) continue;
					}
					else
					{
						if (vertex[n][2] < domain[2])
						{
							if (sqrt((vertex[n][0]-domain[3])*(vertex[n][0]-domain[3]) + (vertex[n][2]-domain[2])*(vertex[n][2]-domain[2])) > outer) continue;
						}
						else if (vertex[n][2] > domain[5])
						{
							if (sqrt((vertex[n][0]-domain[3])*(vertex[n][0]-domain[3]) + (vertex[n][2]-domain[5])*(vertex[n][2]-domain[5])) > outer) continue;
						}
						else if (vertex[n][0]-domain[3] > outer) continue;
					}
				}
				else
				{
					if (vertex[n][1] < domain[1])
					{
						if (vertex[n][2] < domain[2])
						{
							if (sqrt((vertex[n][1]-domain[1])*(vertex[n][1]-domain[1]) + (vertex[n][2]-domain[2])*(vertex[n][2]-domain[2])) > outer) continue;
						}
						else if (vertex[n][2] > domain[5])
						{
							if (sqrt((vertex[n][1]-domain[1])*(vertex[n][1]-domain[1]) + (vertex[n][2]-domain[5])*(vertex[n][2]-domain[5])) > outer) continue;
						}
						else if (domain[1]-vertex[n][1] > outer) continue;
					}
					else if (vertex[n][1] > domain[4])
					{
						if (vertex[n][2] < domain[2])
						{
							if (sqrt((vertex[n][1]-domain[4])*(vertex[n][1]-domain[4]) + (vertex[n][2]-domain[2])*(vertex[n][2]-domain[2])) > outer) continue;
						}
						else if (vertex[n][2] > domain[5])
						{
							if (sqrt((vertex[n][1]-domain[4])*(vertex[n][1]-domain[4]) + (vertex[n][2]-domain[5])*(vertex[n][2]-domain[5])) > outer) continue;
						}
						else if (vertex[n][1]-domain[4] > outer) continue;
					}
					else if (vertex[n][2]-domain[5] > outer || domain[2]-vertex[n][2] > outer) continue;
				}

				if (sqrt((corner[0][0]-vertex[n][0])*(corner[0][0]-vertex[n][0]) + (corner[0][1]-vertex[n][1])*(corner[0][1]-vertex[n][1]) + (corner[0][2]-vertex[n][2])*(corner[0][2]-vertex[n][2])) < inner && sqrt((corner[1][0]-vertex[n][0])*(corner[1][0]-vertex[n][0]) + (corner[1][1]-vertex[n][1])*(corner[1][1]-vertex[n][1]) + (corner[1][2]-vertex[n][2])*(corner[1][2]-vertex[n][2])) < inner && sqrt((corner[2][0]-vertex[n][0])*(corner[2][0]-vertex[n][0]) + (corner[2][1]-vertex[n][1])*(corner[2][1]-vertex[n][1]) + (corner[2][2]-vertex[n][2])*(corner[2][2]-vertex[n][2])) < inner && sqrt((corner[3][0]-vertex[n][0])*(corner[3][0]-vertex[n][0]) + (corner[3][1]-vertex[n][1])*(corner[3][1]-vertex[n][1]) + (corner[3][2]-vertex[n][2])*(corner[3][2]-vertex[n][2])) < inner && sqrt((corner[4][0]-vertex[n][0])*(corner[4][0]-vertex[n][0]) + (corner[4][1]-vertex[n][1])*(corner[4][1]-vertex[n][1]) + (corner[4][2]-vertex[n][2])*(corner[4][2]-vertex[n][2])) < inner && sqrt((corner[5][0]-vertex[n][0])*(corner[5][0]-vertex[n][0]) + (corner[5][1]-vertex[n][1])*(corner[5][1]-vertex[n][1]) + (corner[5][2]-vertex[n][2])*(corner[5][2]-vertex[n][2])) < inner && sqrt((corner[6][0]-vertex[n][0])*(corner[6][0]-vertex[n][0]) + (corner[6][1]-vertex[n][1])*(corner[6][1]-vertex[n][1]) + (corner[6][2]-vertex[n][2])*(corner[6][2]-vertex[n][2])) < inner && sqrt((corner[7][0]-vertex[n][0])*(corner[7][0]-vertex[n][0]) + (corner[7][1]-vertex[n][1])*(corner[7][1]-vertex[n][1]) + (corner[7][2]-vertex[n][2])*(corner[7][2]-vertex[n][2])) < inner) continue; // domain lies within inner sphere

				rdom = 0.5 * sqrt((domain[3]-domain[0])*(domain[3]-domain[0]) + (domain[4]-domain[1])*(domain[4]-domain[1]) + (domain[5]-domain[2])*(domain[5]-domain[2]));
				dist = sqrt((0.5*domain[0]+0.5*domain[3]-vertex[n][0])*(0.5*domain[0]+0.5*domain[3]-vertex[n][0]) + (0.5*domain[1]+0.5*domain[4]-vertex[n][1])*(0.5*domain[1]+0.5*domain[4]-vertex[n][1]) + (0.5*domain[2]+0.5*domain[5]-vertex[n][2])*(0.5*domain[2]+0.5*domain[5]-vertex[n][2]));

				if (dist <= rdom) // vertex lies within domain enclosing sphere
				{
					n++;
					continue;
				}

				if (((0.5*domain[0]+0.5*domain[3]-vertex[n][0])*lightcone.direction[0] + (0.5*domain[1]+0.5*domain[4]-vertex[n][1])*lightcone.direction[1] + (0.5*domain[2]+0.5*domain[5]-vertex[n][2])*lightcone.direction[2]) / dist >= lightcone.opening) // center of domain lies within opening
				{
					n++;
					continue;
				}

				if (dist > outer && acos(((0.5*domain[0]+0.5*domain[3]-vertex[n][0])*lightcone.direction[0] + (0.5*domain[1]+0.5*domain[4]-vertex[n][1])*lightcone.direction[1] + (0.5*domain[2]+0.5*domain[5]-vertex[n][2])*lightcone.direction[2]) / dist) - acos(lightcone.opening) <= acos((outer*outer + dist*dist - rdom*rdom) / (2. * outer * dist))) // enclosing sphere within opening
				{
					n++;
					continue;
				}

				if (dist <= outer && acos(((0.5*domain[0]+0.5*domain[3]-vertex[n][0])*lightcone.direction[0] + (0.5*domain[1]+0.5*domain[4]-vertex[n][1])*lightcone.direction[1] + (0.5*domain[2]+0.5*domain[5]-vertex[n][2])*lightcone.direction[2]) / dist) - acos(lightcone.opening) <= asin(rdom / dist)) // enclosing sphere within opening
				{
					n++;
				}
			}
		}
	}

	return n;
}


//////////////////////////
// hourMinSec
//////////////////////////
// Description:
//   generates formatted output for cpu-time: hh..h:mm:ss.s
//
// Arguments:
//   seconds    number of seconds
//
// Returns:
//   formatted string
//
//////////////////////////

string hourMinSec(double seconds)
{
	string output;
	char ptr[20];
	int h, m, s, f;

	h = (int) floor(seconds / 3600.);
	seconds -= 3600. * h;
	m = (int) floor(seconds / 60.);
	seconds -= 60. * m;
	s = (int) floor(seconds);
	seconds -= s;
	f = (int) floor(10. * seconds);
	sprintf(ptr, "%d:%02d:%02d.%d", h, m, s, f);

	output.reserve(20);
	output.assign(ptr);

	return output;
}

#endif
