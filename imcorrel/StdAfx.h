// stdafx.h : include file for standard system include files,
//  or project specific include files that are used frequently, but
//      are changed infrequently
//

#if !defined(AFX_STDAFX_H__A24A22A3_A2AB_42BA_937B_7C5A1431877C__INCLUDED_)
#define AFX_STDAFX_H__A24A22A3_A2AB_42BA_937B_7C5A1431877C__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#ifdef _DEBUG
#include <cassert>
#define ASSERT(x) assert(x)
#else
#define ASSERT(x)
#define TRACE(x)
#endif

#include <math.h>
#include <memory>

// Matlab headers
#include <mex.h>
#include <matrix.h>

//{{AFX_INSERT_LOCATION}}
// Microsoft Visual C++ will insert additional declarations immediately before the previous line.

#endif // !defined(AFX_STDAFX_H__A24A22A3_A2AB_42BA_937B_7C5A1431877C__INCLUDED_)
