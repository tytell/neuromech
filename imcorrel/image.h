//////////////////////////////////////////////
// class mxRect - Just a simple rectangle class

class mxRect
{
public:
	mwIndex left,top;
	mwIndex right,bottom;

public:
	inline void SetCenter(mwIndex x, mwIndex y, mwSize w, mwSize h);
	inline void Offset(mwIndex dx, mwIndex dy);

	mwSize Width() { return right - left; }
	mwSize Height() { return bottom - top; }
	mwSize Area() { return (right - left)*(bottom - top); }
};

// mxRect inline functions
void mxRect::Offset(mwIndex dx, mwIndex dy)
{
	left += dx;
	right += dx;

	top += dy;
	bottom += dy;
}

void mxRect::SetCenter(mwIndex x, mwIndex y, mwSize w, mwSize h)
{
	left = x - (mwIndex)w/2;
	top = y - (mwIndex)h/2;

	right = left + w;
	bottom = top + h;
}

/////////////////////////////////////////////////////////////////////
// class mxCorrelation - Stores the correlation matrices at each point from
// the mxGrid.  Also keeps track of errors and ongoing calculation

enum { CT_Statistical = 0, CT_Error };		// correlation types
enum ErrorCorrType { NoError = 0, EC_ZeroVariance = 1, NERRCORR };

class mxCorrelation
{
public:
	// NB: all of these arrays are referenced in column major order (Matlab/Fortran
	// default) rather than row major (C default)

	double *x1, *y1;		// positions in the first image
	double *x2, *y2;		// positions in the second image (may be NULL)
	mwSize npt, *dims, ndim;	// number of points in and size of the (x,y) matrices

	mwSize wWidth, wHeight;	// interrogation window size
	mwSize sWidth, sHeight;	// search window size

protected:
	double **phi;				// array of correlation matrices.  Some may be NULL
	double *phi1;				// matrix for ongoing calculation of a correlation matrix
	ErrorCorrType *error;		// array of errors

	static double mxNaN;		// just a value equivalent to Matlab's NaN

public:
	mwSize pr,pc;					// width, height of an individual correlation matrix (i.e. phi1)
	mwSize nphi;					// total size of an individual correlation matrix
	mwIndex ox1, oy1;				// initial displacement of the search region

public:
	mxCorrelation(const mxArray *arr[]);
	~mxCorrelation();

	bool IsDef(mwIndex k);
	mwIndex X1(mwIndex k) { return (mwIndex)x1[k]-1; }
	mwIndex Y1(mwIndex k) { return (mwIndex)y1[k]-1; }
	mwIndex X2(mwIndex k) { if (x2 == NULL) return (mwIndex)x1[k]-1; else return (mwIndex)x2[k]-1; }
	mwIndex Y2(mwIndex k) { if (y2 == NULL) return (mwIndex)y1[k]-1; else return (mwIndex)y2[k]-1; }

	void SavePhi(mwIndex k);		// moves phi1 into the phi array and allocates a new phi1
	
	void SetPhiError(mwIndex k, ErrorCorrType err)			// sets a value in error
		{ ASSERT((k >= 0) && (k < npt)); error[k] = err; }
	
	double &RC(mwIndex r, mwIndex c)	// column major reference into phi1
		{ ASSERT((r >= 0) && (r < pr) && (c >= 0) && (c < pc)); return phi1[c*pr + r]; }

	// convert our mxCorrelation to an mxArray
	mxArray *GetPhiMxArray();
	// convert the error matrix
	mxArray *GetErrorMxArray();

 private:
	// no copying allowed
	mxCorrelation(const mxCorrelation &corr)	{ }
};

/////////////////////////////////////////////////////////////////
// class mxImage - Stores image data and actually performs the correlation.
// Includes helper functions for adding up rows and columns
class mxImage
{
protected:
	double *image;
	mwSize width, height;

public:
	mxImage(const mxArray *arr);

	// Cross correlate this image with I2 at the points in grid
	void Correlate(mxCorrelation &corr, mxImage *I2, int type = CT_Statistical);

protected:
	double *C(mwIndex c)							// column reference
		{ ASSERT((c >= 0) && (c < width)); return &image[c*height]; }
	double *RC(mwIndex r, mwIndex c)					// row and column reference
		{ ASSERT((r >= 0) && (r < height) && (c >= 0) && (c < width)); return &image[c*height + r]; }

	double rectSum(mxRect rect);				// add up a rectangular region
	double colSum(mwIndex x, mwIndex y1, mwIndex y2);		// sum only a column
	double rectVar(mxRect rect, double mean);	// calculate sum of the squared deviation
};

