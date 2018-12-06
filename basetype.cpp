//#include <algorithm>
#include <math.h>
#include <assert.h>
#include <stdlib.h>
#include "basetype.h"


#define MAX(x, y)  ((x) > (y) ? (x) : (y))
#define MIN(x, y)  ((x) < (y) ? (x) : (y))


#ifdef MULTI_THREAD
// atomic add, thread safe, we use arm instructions to implement it
#define XADD(addr, delta) __atomic_fetch_add((_Atomic(int)*)(addr), (delta), 5)
#else
// warning: this version is not thread safe, if you use matrix in multi-thread, is could be error
static inline int XADD(int* addr, int delta)
{ 
	int tmp = *addr; 
	*addr += delta; 
	return tmp; 
}
#endif

static inline size_t alignSize(size_t sz, int n)
{
	//assert((n & (n - 1)) == 0); // n is a power of 2
	return (sz + n-1) & -n;
}

/////////////////////////////CPoint///////////////////////////////////////////////////////////////////////////////////////////////////////////
CPoint::CPoint() : x(0), y(0) 
{
}

CPoint::CPoint(int _x, int _y) : x(_x), y(_y) 
{
}

CPoint::CPoint(const CPoint& pt) : x(pt.x), y(pt.y) 
{
}

CPoint::CPoint(const CSize& sz) : x(sz.width), y(sz.height) 
{
}

CPoint& CPoint::operator = (const CPoint& pt)
{ 
	x = pt.x; 
	y = pt.y; 
	return *this; 
}

static CPoint& operator += (CPoint& a, const CPoint& b)
{
	a.x = a.x + b.x;
	a.y = a.y + b.y;
	return a;
}

static CPoint& operator -= (CPoint& a, const CPoint& b)
{
	a.x = a.x - b.x;
	a.y = a.y - b.y;
	return a;
}

static CPoint& operator *= (CPoint& a, int b)
{
	a.x = a.x*b;
	a.y = a.y*b;
	return a;
}

static CPoint& operator *= (CPoint& a, float b)
{
	a.x = static_cast<int>(a.x*b);
	a.y = static_cast<int>(a.y*b);
	return a;
}

static CPoint& operator *= (CPoint& a, double b)
{
	a.x = static_cast<int>(a.x*b);
	a.y = static_cast<int>(a.y*b);
	return a;
}

static double norm(const CPoint& pt)
{ 
	//return std::sqrt((double)pt.x*pt.x + (double)pt.y*pt.y); 
	return sqrt((double)pt.x*pt.x + (double)pt.y*pt.y); 
}

static bool operator == (const CPoint& a, const CPoint& b)
{ 
	return a.x == b.x && a.y == b.y; 
}

static bool operator != (const CPoint& a, const CPoint& b)
{ 
	return a.x != b.x || a.y != b.y; 
}

static CPoint operator + (const CPoint& a, const CPoint& b)
{ 
	return CPoint( static_cast<int>(a.x + b.x), static_cast<int>(a.y + b.y) ); 
}

static CPoint operator - (const CPoint& a, const CPoint& b)
{ 
	return CPoint( static_cast<int>(a.x - b.x), static_cast<int>(a.y - b.y) ); 
}

static CPoint operator - (const CPoint& a)
{ 
	return CPoint( static_cast<int>(-a.x), static_cast<int>(-a.y) ); 
}

static CPoint operator * (const CPoint& a, int b)
{ 
	return CPoint( static_cast<int>(a.x*b), static_cast<int>(a.y*b) ); 
}

static CPoint operator * (int a, const CPoint& b)
{ 
	return CPoint( static_cast<int>(b.x*a), static_cast<int>(b.y*a) ); 
}

 static CPoint operator * (const CPoint& a, float b)
{ 
	return CPoint( static_cast<int>(a.x*b), static_cast<int>(a.y*b) ); 
 }

 static CPoint operator * (float a, const CPoint& b)
{ 
	return CPoint( static_cast<int>(b.x*a), static_cast<int>(b.y*a) ); 
 }

 static CPoint operator * (const CPoint& a, double b)
{ 
	return CPoint( static_cast<int>(a.x*b), static_cast<int>(a.y*b) ); 
 }

 static CPoint operator * (double a, const CPoint& b)
{ 
	return CPoint( static_cast<int>(b.x*a), static_cast<int>(b.y*a) ); 
 }

 CPoint2f::CPoint2f() : x(0.0), y(0.0)
 {
 }

 CPoint2f::CPoint2f(float _x, float _y) : x(_x), y(_y)
 {
 }

 CPoint2f::CPoint2f(const CPoint2f& pt) : x(pt.x), y(pt.y)
 {
 }

 //CPoint2f::CPoint2f(const CSize& sz) : x(sz.width), y(sz.height)
 //{
 //}

 CPoint2f& CPoint2f::operator = (const CPoint2f& pt)
 {
	 x = pt.x;
	 y = pt.y;
	 return *this;
 }
//////////////////////////////// CSize /////////////////////////////////////////////////////////////////////////////
CSize::CSize(): width(0), height(0) 
{
}

CSize::CSize(int _width, int _height): width(_width), height(_height) 
{
}

CSize::CSize(const CSize& sz): width(sz.width), height(sz.height) 
{
}

CSize::CSize(const CPoint& pt) : width(pt.x), height(pt.y) 
{
}

CSize& CSize::operator = (const CSize& sz)
{ 
	width = sz.width; height = sz.height; return *this; 
}

static CSize operator * (const CSize& a, int b)
{ 
	return CSize(a.width * b, a.height * b); 
}

static CSize operator + (const CSize& a, const CSize& b)
{ 
	return CSize(a.width + b.width, a.height + b.height); 
}

static CSize operator - (const CSize& a, const CSize& b)
{ 
	return CSize(a.width - b.width, a.height - b.height); 
}

int CSize::area() const 
{ 
	return width*height; 
}

static CSize& operator += (CSize& a, const CSize& b)
{ 
	a.width += b.width; a.height += b.height; return a; 
}

static CSize& operator -= (CSize& a, const CSize& b)
{ 
	a.width -= b.width; a.height -= b.height; return a; 
}

static bool operator == (const CSize& a, const CSize& b)
{ 
	return a.width == b.width && a.height == b.height; 
}
static bool operator != (const CSize& a, const CSize& b)
{ 
	return a.width != b.width || a.height != b.height; 
}

//////////////////////////////// CRect ////////////////////////////////
CRect::CRect() : x(0), y(0), width(0), height(0) 
{
}

CRect::CRect(int _x, int _y, int _width, int _height) : x(_x), y(_y), width(_width), height(_height) 
{
}

CRect::CRect(const CRect& r) : x(r.x), y(r.y), width(r.width), height(r.height) 
{
}

CRect::CRect(const CPoint& org, const CSize& sz) :x(org.x), y(org.y), width(sz.width), height(sz.height) 
{
}

CRect::CRect(const CPoint& pt1, const CPoint& pt2)
{
	//x = std::min(pt1.x, pt2.x); 
	x = MIN(pt1.x, pt2.x); 
	//y = std::min(pt1.y, pt2.y);
	y = MIN(pt1.y, pt2.y);
	//width = std::max(pt1.x, pt2.x) - x; 
	width = MAX(pt1.x, pt2.x) - x; 
	//height = std::max(pt1.y, pt2.y) - y;
	height = MAX(pt1.y, pt2.y) - y;
}

CRect& CRect::operator = ( const CRect& r )
{ 
	x = r.x; 
	y = r.y; 
	width = r.width; 
	height = r.height; 
	return *this; 
}

CPoint CRect::tl() const 
{ 
	return CPoint(x,y); 
}

CPoint CRect::br() const 
{ 
	return CPoint(x+width, y+height); 
}

static CRect& operator += ( CRect& a, const CPoint& b )
{ 
	a.x += b.x; 
	a.y += b.y; 
	return a; 
}

static CRect& operator -= ( CRect& a, const CPoint& b )
{ 
	a.x -= b.x; 
	a.y -= b.y; 
	return a; 
}

static CRect& operator += ( CRect& a, const CSize& b )
{ 
	a.width += b.width; 
	a.height += b.height; 
	return a; 
}

static CRect& operator -= ( CRect& a, const CSize& b )
{ 
	a.width -= b.width; 
	a.height -= b.height; 
	return a; 
}

static CRect& operator &= ( CRect& a, const CRect& b )
{
	//int x1 = std::max(a.x, b.x), y1 = std::max(a.y, b.y);
	int x1 = MAX(a.x, b.x), y1 = MAX(a.y, b.y);
	//a.width = std::min(a.x + a.width, b.x + b.width) - x1;
	a.width = MIN(a.x + a.width, b.x + b.width) - x1;
	//a.height = std::min(a.y + a.height, b.y + b.height) - y1;
	a.height = MIN(a.y + a.height, b.y + b.height) - y1;
	a.x = x1; a.y = y1;
	if( a.width <= 0 || a.height <= 0 )
		a = CRect();
	return a;
}

static CRect& operator |= ( CRect& a, const CRect& b )
{
	//int x1 = std::min(a.x, b.x), y1 = std::min(a.y, b.y);
	int x1 = MIN(a.x, b.x), y1 = MIN(a.y, b.y);
	//a.width = std::max(a.x + a.width, b.x + b.width) - x1;
	a.width = MAX(a.x + a.width, b.x + b.width) - x1;
	//a.height = std::max(a.y + a.height, b.y + b.height) - y1;
	a.height = MAX(a.y + a.height, b.y + b.height) - y1;
	a.x = x1; a.y = y1;
	return a;
}

CSize CRect::size() const 
{ 
	return CSize(width, height); 
}

int CRect::area() const 
{ 
	return width*height; 
}

bool CRect::contains(const CPoint& pt) const
{ 
	return x <= pt.x && pt.x < x + width && y <= pt.y && pt.y < y + height; 
}

static bool operator == (const CRect& a, const CRect& b)
{
	return a.x == b.x && a.y == b.y && a.width == b.width && a.height == b.height;
}

static bool operator != (const CRect& a, const CRect& b)
{
	return a.x != b.x || a.y != b.y || a.width != b.width || a.height != b.height;
}

static CRect operator + (const CRect& a, const CPoint& b)
{
	return CRect( a.x + b.x, a.y + b.y, a.width, a.height );
}

static CRect operator - (const CRect& a, const CPoint& b)
{
	return CRect( a.x - b.x, a.y - b.y, a.width, a.height );
}

static CRect operator + (const CRect& a, const CSize& b)
{
	return CRect( a.x, a.y, a.width + b.width, a.height + b.height );
}

static CRect operator & (const CRect& a, const CRect& b)
{
	CRect c = a;
	return c &= b;
}

static CRect operator | (const CRect& a, const CRect& b)
{
	CRect c = a;
	return c |= b;
}

bool CPoint::inside( const CRect& r ) const
{
	return r.contains(*this);
}

CRect2f::CRect2f() : x(0.0), y(0.0), width(0.0), height(0.0)
{
}

CRect2f::CRect2f(float _x, float _y, float _width, float _height) : x(_x), y(_y), width(_width), height(_height)
{
}

CRect2f::CRect2f(const CRect2f& r) : x(r.x), y(r.y), width(r.width), height(r.height)
{
}

CRect2f& CRect2f::operator = (const CRect2f& r)
{
	x = r.x;
	y = r.y;
	width = r.width;
	height = r.height;
	return *this;
}
////////////////////////////////////////////////////////////////////////////////////////////////
// CMat is the matrix class, support int and float, multi-channel, multi-dimensional

// _dims : dimension
// _sz: size of each dimension
// _steps: step of each dimension
static inline void setSize( CMat& m, int _dims, const int* _sz,
	const size_t* _steps, bool autoSteps=false )
{
	//assert( 0 <= _dims && _dims <= MAT_MAX_DIM );
	// if the matrix is not empty, free it first
	if ( m.dims != _dims )
	{
		if ( m.step.p != m.step.buf )
		{
			free(m.step.p);
			m.step.p = m.step.buf;
			m.size.p = &m.rows;
		}

		if ( _dims > 2 )
		{
			m.step.p = (size_t*)malloc(_dims*sizeof(m.step.p[0]) + (_dims+1)*sizeof(m.size.p[0]));
			m.size.p = (int*)(m.step.p + _dims) + 1;
			m.size.p[-1] = _dims; // dim, s1, s2,..., sn
			m.rows = m.cols = -1;
		}
	}

	m.dims = _dims;
	if ( !_sz )
		return;

	size_t esz = MAT_ELEM_SIZE(m.flags), total = esz;
	int i;
	for ( i = _dims - 1; i >= 0; i-- )
	{
		int s = _sz[i];
		//assert( s >= 0 );
		m.size.p[i] = s; // set size for each dimension

		if( _steps )
			m.step.p[i] = i < _dims-1 ? _steps[i] : esz;
		else if( autoSteps )
		{
			m.step.p[i] = total;
			//int64 total1 = (int64)total*s;
			size_t total1 = total*s; // warning: if total*s is too big, it must be error
			//assert( total1 == (size_t)total1 );
			total = (size_t)total1;
		}
	}

	if( _dims == 1 )
	{
		m.dims = 2;
		m.cols = 1;
		m.step[1] = esz;
	}
}

// initialize matrix as empty
void CMat::initEmpty()
{
    flags = MAT_MAGIC_VAL;
    dims = rows = cols = 0;
    data = datastart = dataend = datalimit = 0;
    refcount = 0;
}

// ceate a empty matrix
CMat::CMat() : size(&rows)
{
    initEmpty();
}

// ceate a 2d matrix
CMat::CMat(int _rows, int _cols, int _type) : size(&rows)
{
    initEmpty();
    create(_rows, _cols, _type);
}

// create a 2d matrix
CMat::CMat(CSize _sz, int _type) : size(&rows)
{
    initEmpty();
    create( _sz.height, _sz.width, _type );
}

// copy constructer, data is shared between to matrix
CMat::CMat(const CMat& m)
    : flags(m.flags), dims(m.dims), rows(m.rows), cols(m.cols), data(m.data),
    refcount(m.refcount), datastart(m.datastart), dataend(m.dataend),
    datalimit(m.datalimit), size(&rows)
{
    if ( refcount ) // increment reference counter
        XADD(refcount, 1);
    if ( m.dims <= 2 )
    {
        step[0] = m.step[0]; step[1] = m.step[1];
    }
    else
    {
        dims = 0;
        copySize(m);
    }
}

CMat::CMat(int _rows, int _cols, int _type, void* _data, size_t _step)
    : flags(MAT_MAGIC_VAL + (_type & MAT_TYPE_MASK)), dims(2), rows(_rows), cols(_cols),
    data((uchar*)_data), refcount(0), datastart((uchar*)_data), dataend(0),
    datalimit(0), size(&rows)
{
    size_t esz = MAT_ELEM_SIZE(_type), minstep = cols*esz;
    if( _step == CMAT_AUTO_STEP )
    {
        _step = minstep;
        flags |= MAT_CONTINUOUS_FLAG;
    }
    else
    {
        if ( rows == 1 ) 
			_step = minstep;
        //assert( _step >= minstep );
        flags |= _step == minstep ? MAT_CONTINUOUS_FLAG : 0;
    }
    step[0] = _step; step[1] = esz;
    datalimit = datastart + _step*rows;
    dataend = datalimit - _step + minstep;
}

CMat::CMat(CSize _sz, int _type, void* _data, size_t _step)
    : flags(MAT_MAGIC_VAL + (_type & MAT_TYPE_MASK)), dims(2), rows(_sz.height), cols(_sz.width),
    data((uchar*)_data), refcount(0), datastart((uchar*)_data), dataend(0),
    datalimit(0), size(&rows)
{
    size_t esz = MAT_ELEM_SIZE(_type), minstep = cols*esz;
    if ( _step == CMAT_AUTO_STEP )
    {
        _step = minstep;
        flags |= MAT_CONTINUOUS_FLAG;
    }
    else
    {
        if ( rows == 1 ) 
			_step = minstep;
        //assert( _step >= minstep );
        flags |= _step == minstep ? MAT_CONTINUOUS_FLAG : 0;
    }
    step[0] = _step; step[1] = esz;
    datalimit = datastart + _step*rows;
    dataend = datalimit - _step + minstep;
}

// release data if the reference counter is 0
CMat::~CMat()
{
    release();
    if ( step.p != step.buf )
	    free(step.p);
}

CMat& CMat::operator = (const CMat& m)
{
    if ( this != &m )
    {
		if ( m.refcount )
			XADD(m.refcount, 1);
        release();
        flags = m.flags;
        if ( dims <= 2 && m.dims <= 2 )
        {
            dims = m.dims;
            rows = m.rows;
            cols = m.cols;
            step[0] = m.step[0];
            step[1] = m.step[1];
        }
        else
            copySize(m);
        data = m.data;
        datastart = m.datastart;
        dataend = m.dataend;
        datalimit = m.datalimit;
        refcount = m.refcount;
    }
    return *this;
}

// create a 2d matrix
void CMat::create(int _rows, int _cols, int _type)
{
    _type &= MAT_TYPE_MASK;
    if ( dims <= 2 && rows == _rows && cols == _cols && type() == _type && data )
        return;
    int sz[] = {_rows, _cols};
    create(2, sz, _type);
}

static void finalizeHdr(CMat& m)
{
	//updateContinuityFlag(m);
	int d = m.dims;
	if ( d > 2 )
		m.rows = m.cols = -1;
	if ( m.data )
	{
		m.datalimit = m.datastart + m.size[0]*m.step[0];
		if ( m.size[0] > 0 )
		{
			m.dataend = m.data + m.size[d-1]*m.step[d-1];
			for ( int i = 0; i < d-1; i++ )
				m.dataend += (m.size[i] - 1)*m.step[i];
		}
		else
			m.dataend = m.datalimit;
	}
	else
		m.dataend = m.datalimit = 0;
}

// create a n-d matrix
void CMat::create(int d, const int* _sizes, int _type)
{
	int i;
	//assert(0 <= d && d <= MAT_MAX_DIM && _sizes);
	//_type = MAT_MAT_TYPE(_type);

	if ( data && (d == dims || (d == 1 && dims <= 2)) && _type == type() )
	{
		if ( d == 2 && rows == _sizes[0] && cols == _sizes[1] )
			return;
		for ( i = 0; i < d; i++ )
			if ( size[i] != _sizes[i] )
				break;
		if ( i == d && (d > 1 || size[1] == 1))
			return;
	}

	release();
	if ( d == 0 )
		return;
	flags = (_type & MAT_MAT_TYPE_MASK) | MAT_MAGIC_VAL;
	setSize(*this, d, _sizes, 0, true);

	if ( total() > 0 )
	{
		size_t totalsize = alignSize(step.p[0]*size.p[0], (int)sizeof(*refcount));
		data = datastart = (uchar*)malloc(totalsize + (int)sizeof(*refcount));
		refcount = (int*)(data + totalsize);
		*refcount = 1;
	}

	finalizeHdr(*this);
}

// create a 2d matrx
void CMat::create(CSize _sz, int _type)
{
    create(_sz.height, _sz.width, _type);
}

void CMat::addref()
{ 
	if ( refcount ) 
		XADD(refcount, 1);
}

void CMat::release()
{
	if ( refcount && XADD(refcount, -1) == 1 )
        deallocate();
    data = datastart = dataend = datalimit = 0;
    size.p[0] = 0;
    refcount = 0;
}

void CMat::deallocate()
{
	//assert(refcount != 0);
	free(datastart);
}

void CMat::copySize(const CMat& m)
{
	setSize(*this, m.dims, 0, 0);
	for ( int i = 0; i < dims; i++ )
	{
		size[i] = m.size[i];
		step[i] = m.step[i];
	}
}

int CMat::type() const 
{ 
	return MAT_MAT_TYPE(flags); 
}

int CMat::depth() const 
{ 
	return MAT_MAT_DEPTH(flags); 
}

int CMat::channels() const 
{ 
	return MAT_MAT_CN(flags); 
}

//size_t CMat::step1(int i) const 
//{ 
//	return step.p[i]/elemSize1(); 
//}

bool CMat::empty() const 
{ 
	return data == 0 || total() == 0; 
}

size_t CMat::total() const
{
    if ( dims <= 2 )
        return (size_t)rows*cols;
    size_t p = 1;
    for ( int i = 0; i < dims; i++ )
        p *= size[i];
    return p;
}

uchar* CMat::ptr(int y)
{
    //assert( y == 0 || (data && dims >= 1 && (unsigned)y < (unsigned)size.p[0]) );
    return data + step.p[0]*y;
}

const uchar* CMat::ptr(int y) const
{
    //assert( y == 0 || (data && dims >= 1 && (unsigned)y < (unsigned)size.p[0]) );
    return data + step.p[0]*y;
}

uchar* CMat::ptr(int i0, int i1)
{
    //assert( dims >= 2 && data &&
    //              (unsigned)i0 < (unsigned)size.p[0] &&
    //              (unsigned)i1 < (unsigned)size.p[1] );
    return data + i0*step.p[0] + i1*step.p[1];
}

const uchar* CMat::ptr(int i0, int i1) const
{
    //assert( dims >= 2 && data &&
    //             (unsigned)i0 < (unsigned)size.p[0] &&
    //             (unsigned)i1 < (unsigned)size.p[1] );
    return data + i0*step.p[0] + i1*step.p[1];
}

uchar* CMat::ptr(int i0, int i1, int i2)
{
	//assert( dims >= 3 && data &&
	//	(unsigned)i0 < (unsigned)size.p[0] &&
	//	(unsigned)i1 < (unsigned)size.p[1] &&
	//	(unsigned)i2 < (unsigned)size.p[2] );
    return data + i0*step.p[0] + i1*step.p[1] + i2*step.p[2];
}

const uchar* CMat::ptr(int i0, int i1, int i2) const
{
    //assert( dims >= 3 && data &&
    //              (unsigned)i0 < (unsigned)size.p[0] &&
    //              (unsigned)i1 < (unsigned)size.p[1] &&
    //              (unsigned)i2 < (unsigned)size.p[2] );
    return data + i0*step.p[0] + i1*step.p[1] + i2*step.p[2];
}

uchar* CMat::ptr(const int* idx)
{
    int i, d = dims;
    uchar* p = data;
    //assert( d >= 1 && p );
    for ( i = 0; i < d; i++ )
    {
        //assert( (unsigned)idx[i] < (unsigned)size.p[i] );
        p += idx[i]*step.p[i];
    }
    return p;
}

const uchar* CMat::ptr(const int* idx) const
{
    int i, d = dims;
    uchar* p = data;
    //assert( d >= 1 && p );
    for( i = 0; i < d; i++ )
    {
        //assert( (unsigned)idx[i] < (unsigned)size.p[i] );
        p += idx[i]*step.p[i];
    }
    return p;
}

//2017.4.16 jimmypan added
CMat CMat::mul(CMat &m){
	assert(cols == m.cols && rows == m.rows);
	CMat res = CMat(rows, cols, MAT_32F);
	for (int i = 0; i < rows; i++)
	{
		float* pres = (float*)res.ptr(i);
		float* pth = (float*)ptr(i);
		float* pm = (float*)m.ptr(i);
		for (int j = 0; j < cols; j++)
		{
			pres[j] = pth[j] * pm[j];
		}
	}
	return res;
}
CMat operator * (CMat m,float fval){
	CMat res = CMat(m.rows, m.cols, m.type());
	int ch = m.channels();
	for (int i = 0; i < m.rows; i++)
	{
		float* pres = (float*)res.ptr(i);
		float* pth = (float*)m.ptr(i);
		for (int j = 0; j < m.cols; j++)
		{
			for (int k = 0; k < ch; k++)
			{
				int idx = j*ch + k;
				pres[idx] = pth[idx] * fval;
			}			
		}
	}
	return res;
}
CMat operator + (CMat m,float fval){
	CMat res = CMat(m.rows, m.cols, m.type());
	int ch = m.channels();
	for (int i = 0; i < m.rows; i++)
	{
		float* pres = (float*)res.ptr(i);
		float* pth = (float*)m.ptr(i);
		for (int j = 0; j < m.cols; j++)
		{
			int idx = j*ch;
			pres[idx] = pth[idx] + fval;
			for (int k = 1; k < ch; k++)
			{
				pres[idx + k] = pth[idx + k];
			}
		}
	}
	return res;
}
CMat operator + (CMat m, CMat n){
	assert(n.cols == m.cols && n.rows == m.rows && n.type() == m.type());
	CMat res = CMat(m.rows, m.cols, m.type());
	int ch = m.channels();
	for (int i = 0; i < m.rows; i++)
	{
		float* pres = (float*)res.ptr(i);
		float* pth = (float*)n.ptr(i);
		float* pm = (float*)m.ptr(i);
		for (int j = 0; j < m.cols; j++)
		{
			for (int k = 0; k < ch; k++)
			{
				int idx = j*ch + k;
				pres[idx] = pth[idx] + pm[idx];
			}			
		}
	}
	return res;
}

//////////////////////////////////////////////////////////////////////////////
// MSize
CMat::MSize::MSize(int* _p) : p(_p) 
{
}

CSize CMat::MSize::operator()() const
{
    //assert(p[-1] <= 2);
    return CSize(p[1], p[0]);
}

const int& CMat::MSize::operator[](int i) const 
{ 
	return p[i]; 
}

int& CMat::MSize::operator[](int i) 
{ 
	return p[i]; 
}

CMat::MSize::operator const int*() const 
{ 
	return p; 
}

bool CMat::MSize::operator == (const MSize& sz) const
{
    int d = p[-1], dsz = sz.p[-1];
    if ( d != dsz )
        return false;
    if ( d == 2 )
        return p[0] == sz.p[0] && p[1] == sz.p[1];

    for ( int i = 0; i < d; i++ )
        if ( p[i] != sz.p[i] )
            return false;
    return true;
}

bool CMat::MSize::operator != (const MSize& sz) const
{
    return !(*this == sz);
}

///////////////////////////////////////////////////////////////////////////
// MStep
CMat::MStep::MStep() 
{ 
	p = buf; 
	p[0] = p[1] = 0; 
}

CMat::MStep::MStep(size_t s) 
{ 
	p = buf; 
	p[0] = s; 
	p[1] = 0; 
}

const size_t& CMat::MStep::operator[](int i) const 
{ 
	return p[i]; 
}

size_t& CMat::MStep::operator[](int i) 
{ 
	return p[i]; 
}

CMat::MStep::operator size_t() const
{
    //assert( p == buf );
    return buf[0];
}

CMat::MStep& CMat::MStep::operator = (size_t s)
{
    //assert( p == buf );
    buf[0] = s;
    return *this;
}
