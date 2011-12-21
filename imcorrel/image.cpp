#include "stdafx.h"
#include "imcorrel.h"
#include "image.h"

long longnan[2] = { 0x00000000, 0xfff80000 };
double *dbnan = (double *)&longnan;
double mxCorrelation::mxNaN = *dbnan;

/////////////////////////////////////////////////////////////
// mxImage::Correlate - Cross-correlate this image with I2 at the points in
// grid, using the type of correlation specified in type.
//
// TODO: Add error correlation.  Currently only statistical correlation is supported
void mxImage::Correlate(mxCorrelation &corr, mxImage *I2, int type)
{
	mwIndex k;
	mwIndex ox, oy;
	mwIndex wc,sc, r;
	double *wCol, *sCol;
	mwIndex x1,y1, x2,y2;
	double wMean,wVar;
	double sSum,sMean;
	double sVar;
	double corrVal;

	mxRect win, search;

	for (k = 0; k < corr.npt; k++) {
		// use the try/catch block here to exit gracefully from inner loops and end
		// up back here, without using a goto or excessive flags
		// Note that we don't handle our general error type ErrorType here, so severe
		// errors get passed up to the global try/catch block
		try {
			if (!corr.IsDef(k))
				throw NoError;			// skip this point, but no error

			// get the window position
			x1 = (mwIndex)corr.X1(k);
			y1 = (mwIndex)corr.Y1(k);
			win.SetCenter(x1,y1, corr.wWidth, corr.wHeight);
			if ((win.left < 0) || (win.right >= width) || (win.top < 0) || (win.bottom >= height))
				throw E_GridOutOfBounds;		// this is a real error

			// store its mean and variance for later
			wMean = rectSum(win) / win.Area();
			wVar = rectVar(win, wMean);
			if (wVar == 0)
				throw EC_ZeroVariance;

			// get the search position -- may be displaced from the window position
			// if we're doing superresolution PIV
			x2 = (mwIndex)corr.X2(k);
			y2 = (mwIndex)corr.Y2(k);

			mwIndex minx,miny;
			minx = x2 - corr.sWidth/2;
			miny = y2 - corr.sHeight/2;
			if ((minx < 0) || (minx+corr.sWidth >= width) || 
					(miny < 0) || (miny+corr.sHeight >= height))
				throw E_GridOutOfBounds;

			for (oy = 0; oy < corr.pr; oy++) {
				// NB: the size of the search rect is (wWidth,wHeight), and then
				// we offset it around until the total area searched is (sWidth,
				// sHeight)
				search.SetCenter(x2+corr.ox1, y2+oy+corr.oy1, corr.wWidth, corr.wHeight);

				// get the initial sum over the region
				sSum = I2->rectSum(search);
				sMean = sSum/search.Area();

				for (ox = 0; ox < corr.pc; ox++) {
					sVar = I2->rectVar(search, sMean);
					if (sVar == 0)
						throw EC_ZeroVariance;

					// calculate the unnormalized correlation
					corrVal = 0;
					for (wc = win.left, sc = search.left; wc < win.right; wc++, sc++) {
						wCol = RC(win.top, wc);
						sCol = I2->RC(search.top, sc);
						for (r = 0; r < corr.wHeight; r++) {
							corrVal += (*wCol++ - wMean)*(*sCol++ - sMean);
						}
					}

					// normalize by variance
					corrVal /= sqrt(wVar*sVar);
					// and save this calculated value in temporary matrix (phi1) in corr
					corr.RC(oy,ox) = corrVal;

					// check to make sure we're reasonable
					ASSERT((corrVal >= -1) && (corrVal <= 1));

					// rather than adding up the whole mean each time, we realize that
					// by just moving one column over, the majority of the sum is the same.
					// just subtract the left column and add in a new right one, and we have
					// the new sum
					double firstSum, lastSum;

					firstSum = I2->colSum(search.left, search.top,search.bottom);
					lastSum = I2->colSum(search.right, search.top,search.bottom);
					sSum += lastSum - firstSum;
					sMean = sSum/search.Area();

					// offset the search region
					search.Offset(1,0);
					ASSERT(fabs(sSum - I2->rectSum(search)) < 0.0001);
				}
			}

			// save the temporary correlation matrix (phi1) in the right position in the
			// global matrix
			corr.SavePhi(k);
		}
		catch (ErrorCorrType err) {
			if (err != NoError)
				corr.SetPhiError(k, err);
			// continue with the loop
		}
	}
}

void mxCorrelation::SavePhi(mwIndex k)
// Save the current phi1 matrix in the global phi matrix at position k and
// allocate a new phi1 matrix
{
	phi[k] = phi1;

	// in case the new fails, we need to know that phi1 is in the phi matrix, so
	// we set it to NULL
	phi1 = NULL;
	phi1 = new double[pr*pc];
}

bool mxCorrelation::IsDef(mwIndex k)
// check to see if a grid point is non-NaN
{
	ASSERT((k >= 0) && (k < npt));
	if ((x2 != NULL) && (!mxIsFinite(x2[k]) || !mxIsFinite(y2[k])))
		return false;
	if (!mxIsFinite(x1[k]) || !mxIsFinite(y1[k]))
		return false;
	return true;
}

// mxRect operations
double mxImage::rectSum(mxRect rect)
{
	ASSERT((rect.left >= 0) && (rect.right < width) &&
		(rect.top >= 0) && (rect.bottom < height));

	mwIndex r,c;
	double *col;
	double sum = 0;

	for (c = rect.left; c < rect.right; c++) {
		col = RC(rect.top, c);

		for (r = rect.top; r < rect.bottom; r++) {
			sum += *col++;
		}
	}
	
	return sum;
}

double mxImage::colSum(mwIndex x, mwIndex y1, mwIndex y2)
{
	ASSERT((x >= 0) && (x < width) &&
		(y1 >= 0) && (y2 < height));
	
	mwIndex r;
	double *col;
	double sum = 0;

	col = RC(y1, x);
	for (r = y1; r < y2; r++)
		sum += *col++;

	return sum;
}

double mxImage::rectVar(mxRect rect, double mean)
{
	ASSERT((rect.left >= 0) && (rect.right < width) &&
		(rect.top >= 0) && (rect.bottom < height));

	mwIndex r,c;
	double *col;
	double mb, var = 0;

	for (c = rect.left; c < rect.right; c++) {
		col = RC(rect.top, c);

		for (r = rect.top; r < rect.bottom; r++) {
			mb = *col++ - mean;
			var = var + mb*mb;
		}
	}

	return var;
}		

// Constructors
mxImage::mxImage(const mxArray *arr)
{
	image = mxGetPr(arr);
	width = mxGetN(arr);
	height = mxGetM(arr);
}

mxCorrelation::mxCorrelation(const mxArray *arr[])
{
	// we set these to NULL initially so that if we fail on one of the allocations
	// we'll know which arrays are valid and which aren't
	phi = NULL;
	phi1 = NULL;
	error = NULL;
	dims = NULL;

	// get the grid

	const mwSize *mxDims;
	mwIndex i;

	// get the positions for the first image
	// NB: these are all referenced with the first position being 1, not 0.  We use
	// X1,X2, Y1, and Y2 to get everything in the proper form
	x1 = mxGetPr(arr[0]);
	y1 = mxGetPr(arr[1]);

	// save their size.  We'll reference them just by (Matlab-style) index throughout
	// imcorrel, until we have to return them
	ndim = mxGetNumberOfDimensions(arr[0]);
	dims = new mwSize[ndim];

	mxDims = mxGetDimensions(arr[0]);
	memcpy(dims, mxDims, ndim*sizeof(*dims));

	// if the last dimension is 2, it means we have different positions for the first
	// and second images.  Set (x2,y2) to the second positions if they exist
	for (i = 0, npt = 1; i < ndim; i++)
		npt *= dims[i];

	if (dims[ndim-1] == 2) {
		npt = npt/2;

		x2 = &x1[npt];
		y2 = &y1[npt];
	}
	else {
		x2 = NULL;		// or NULL if they don't
		y2 = NULL;
	}

	// save the window sizes
	double *ptr;

	ptr = mxGetPr(arr[2]);
	wWidth = (mwSize)ptr[0];
	wHeight = (mwSize)ptr[1];

	ptr = mxGetPr(arr[3]);
	sWidth = (mwSize)ptr[0];
	sHeight = (mwSize)ptr[1];
	
	pr = sHeight - wHeight;
	pc = sWidth - wWidth;
	nphi = pr*pc;

	ox1 = -pc/2;
	oy1 = -pr/2;

	// now allocate things

	// allocate
	phi = new double * [npt];
	error = new ErrorCorrType[npt];

	// and set to 0/NULL
	memset(phi, 0, npt*sizeof(double *));
	memset(error, 0, npt*sizeof(ErrorCorrType));

	// temporary storage
	phi1 = new double[nphi];
}

mxCorrelation::~mxCorrelation()
{
	TRACE("mxCorrelation::~mxCorrelation\n");

	if (phi != NULL) {
		mwIndex i;
		for (i = 0; i < npt; i++) {
			if (phi[i] != NULL)
				delete phi[i];
		}
		delete phi;
	}
	if (phi1 != NULL)
		delete phi1;
	if (error != NULL)
		delete error;
	if (dims != NULL)
		delete dims;
}
