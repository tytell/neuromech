// imcorrel.cpp : Defines the entry point for the DLL application.
//

#include "stdafx.h"
#include "imcorrel.h"
#include "image.h"

// #define WIN32_LEAN_AND_MEAN		// Exclude rarely-used stuff from Windows headers
// #include <windows.h>

void CheckParameters(int nrhs, const mxArray *prhs[]);
extern char *ErrorText[NERR];

extern void _main();

int memory_error_handler( size_t )
{
	throw E_OutOfMemory;
}

// imcorrel(I1,I2, x,y, win,search, ...)

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, 
     const mxArray *prhs[])
{
    CheckParameters(nrhs, prhs);
    try {
	mxImage I1(prhs[0]);
	mexPrintf("Here!\n");

	mxImage I2(prhs[1]);

	mxCorrelation phi(&prhs[2]);

	// run the correlation
	I1.Correlate(phi, &I2);

	// return the parameters
	if (nlhs > 0)
	    plhs[0] = phi.GetPhiMxArray();
	if (nlhs == 2)
	    plhs[1] = phi.GetErrorMxArray();
    }
    catch (ErrorType err) {
	mexErrMsgTxt(ErrorText[err]);
    }
    // we have to include this handler under debug code, because we include MFC just
    // enough to use ASSERT and TRACE.  Unfortunately that also overrides our definition
    // of new so that it throws a CException class, not an ErrorType as we want it to.
    // In the Release version, we don't have any MFC, so it doesn't matter
    // #ifdef _DEBUG
    // 	catch (CException *) {
    // 		mexErrMsgTxt("MFC exception");
    // 	}
    // #endif
    catch (...) {
	mexErrMsgTxt("Unknown exception.");
    }
}

/////////////////////////////////////////////
// CheckParameters and various support functions 

inline bool mxIsRealDouble(const mxArray *arr)
{
	return (mxIsDouble(arr) && !mxIsComplex(arr));
}

inline bool mxIsSizeEqual(const mxArray *arr1, const mxArray *arr2)
{
	return (mxGetM(arr1) == mxGetM(arr2)) && (mxGetN(arr1) == mxGetN(arr2));
}

inline bool mxIsSizeEqual(const mxArray *arr, mwSize m, mwSize n)
{
	return (mxGetM(arr) == m) && (mxGetN(arr) == n);
}

void CheckParameters(int nrhs, const mxArray *prhs[])
{
	if (nrhs < 6)
		mexErrMsgTxt("You must supply at least 6 arguments");

	if (!mxIsRealDouble(prhs[0]) || !mxIsRealDouble(prhs[1]))
		mexErrMsgTxt("The images must both be real doubles (values between 0 and 1).");
	if (!mxIsSizeEqual(prhs[0], prhs[1]))
		mexErrMsgTxt("The images must be the same size.");

	if (!mxIsRealDouble(prhs[2]) || !mxIsRealDouble(prhs[3]))
		mexErrMsgTxt("The interrogation window positions must be real values.");

	// the x and y positions need to be matrices with size [a,b,...,2], where they
	// can have any number of dimensions, as long as the last dimension has two elements
	// The first element is the position in the first image, and the second is the position
	// in the second image.  This allows for displacement of the regions between images,
	// which is important in superresolution techniques
	int nd,i;
	const mwSize *dims1, *dims2;

	nd = mxGetNumberOfDimensions(prhs[2]);
	if (mxGetNumberOfDimensions(prhs[3]) != nd)
		mexErrMsgTxt("The interrogation window x and y positions must be the same size.");

	dims1 = mxGetDimensions(prhs[2]);
	dims2 = mxGetDimensions(prhs[3]);
	for (i = 0; i < nd; i++) {
		if (dims1[i] != dims2[i])
			mexErrMsgTxt("The interrogation window x and y positions must be the same size.");
	}

	if ((mxGetNumberOfElements(prhs[4]) != 2) || (mxGetNumberOfElements(prhs[5]) != 2))
		mexErrMsgTxt("The interrogation window sizes must be two-element vectors.");
}

///////////////////////////////////////////////
// Functions for returning the output
//         mxCorrelation::GetPhiMxArray()
//         mxCorrelation::GetErrorMxArray()

mxArray *mxCorrelation::GetPhiMxArray()
{
	mxArray *mxCorrel;
	double *corr;
	mwSize *phidims;
	mwSize nd;
	mwSize i;
	mwSize k, kphi, bytelen;

	nd = ndim + 2;
	phidims = new mwSize [nd];
	phidims[0] = pr;
	phidims[1] = pc;
	for (i = 0; i < nd-2; i++)
		phidims[i+2] = dims[i];

	// ignore the last dimension if it's 1 (i.e. a column vector) or 2 (when they've
	// specified different postions for each image)
	if ((phidims[nd-1] == 1) || (phidims[nd-1] == 2))
		nd = nd-1;

	/* int phidims[3];
	phidims[0] = pr;
	phidims[1] = pc;
	phidims[2] = npt;
	nd = 3; */

	bytelen = nphi * sizeof(double);

	mxCorrel = mxCreateNumericArray(nd,phidims, mxDOUBLE_CLASS, mxREAL);
	corr = mxGetPr(mxCorrel);
	for (k = 0, kphi = 0; k < npt; k++) {
		if (phi[k] != NULL)
			memcpy(&corr[kphi], phi[k], bytelen);
		else {
			double *c = &corr[kphi];
			mwSize j;

			for (j = 0; j < nphi; j++)
				*c++ = mxNaN;	
		}
		kphi += nphi;
	}

	delete phidims;

	return mxCorrel;
}

mxArray *mxCorrelation::GetErrorMxArray()
{
	mxArray *mxError;
	double *matError;
	ErrorCorrType *ecError;
	mwIndex k;
	mwSize *errdims;

	errdims = new mwSize[ndim];
	memcpy(errdims, dims, ndim*sizeof(*errdims));
	if (errdims[ndim-1] == 2)
		errdims[ndim-1] = 1;

	mxError = mxCreateNumericArray(ndim,errdims, mxDOUBLE_CLASS, mxREAL);

	// mxError = mxCreateDoubleMatrix(npt,1, mxREAL);
	matError = mxGetPr(mxError);
	ecError = error;

	// NB: We can't use memcpy, because error is an integer array and matError is a
	// double array
	for (k = 0; k < npt; k++)
		*matError++ = (double)*ecError++;

	delete errdims;

	return mxError;
}
