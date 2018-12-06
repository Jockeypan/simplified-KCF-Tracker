#ifndef __BASE_TYPE_H__
#define __BASE_TYPE_H__

#include <stddef.h>
#include <stdio.h>

//////////////////////////////////////////////////////////////////////////////////////
// constant for CMat
#define  MAT_MAGIC_VAL   0x42FF0000
#define  MAT_TYPE_MASK  0x00000FFF

#define MAT_MAX_DIM            32

#define MAT_CN_MAX     512
#define MAT_CN_SHIFT   3
#define MAT_DEPTH_MAX  (1 << MAT_CN_SHIFT)

#define MAT_8U   0
#define MAT_8S   1
#define MAT_16U  2
#define MAT_16S  3
#define MAT_32S  4
#define MAT_32F  5
#define MAT_64F  6
#define MAT_USRTYPE1 7

#define MAT_MAT_DEPTH_MASK       (MAT_DEPTH_MAX - 1)
#define MAT_MAT_DEPTH(flags)     ((flags) & MAT_MAT_DEPTH_MASK)

#define MAT_MAKETYPE(depth,cn) (MAT_MAT_DEPTH(depth) + (((cn)-1) << MAT_CN_SHIFT))
#define MAT_MAKE_TYPE MAT_MAKETYPE

#define MAT_8UC1 MAT_MAKETYPE(MAT_8U,1)
#define MAT_8UC2 MAT_MAKETYPE(MAT_8U,2)
#define MAT_8UC3 MAT_MAKETYPE(MAT_8U,3)
#define MAT_8UC4 MAT_MAKETYPE(MAT_8U,4)
#define MAT_8UC(n) MAT_MAKETYPE(MAT_8U,(n))

#define MAT_8SC1 MAT_MAKETYPE(MAT_8S,1)
#define MAT_8SC2 MAT_MAKETYPE(MAT_8S,2)
#define MAT_8SC3 MAT_MAKETYPE(MAT_8S,3)
#define MAT_8SC4 MAT_MAKETYPE(MAT_8S,4)
#define MAT_8SC(n) MAT_MAKETYPE(MAT_8S,(n))

#define MAT_16UC1 MAT_MAKETYPE(MAT_16U,1)
#define MAT_16UC2 MAT_MAKETYPE(MAT_16U,2)
#define MAT_16UC3 MAT_MAKETYPE(MAT_16U,3)
#define MAT_16UC4 MAT_MAKETYPE(MAT_16U,4)
#define MAT_16UC(n) MAT_MAKETYPE(MAT_16U,(n))

#define MAT_16SC1 MAT_MAKETYPE(MAT_16S,1)
#define MAT_16SC2 MAT_MAKETYPE(MAT_16S,2)
#define MAT_16SC3 MAT_MAKETYPE(MAT_16S,3)
#define MAT_16SC4 MAT_MAKETYPE(MAT_16S,4)
#define MAT_16SC(n) MAT_MAKETYPE(MAT_16S,(n))

#define MAT_32SC1 MAT_MAKETYPE(MAT_32S,1)
#define MAT_32SC2 MAT_MAKETYPE(MAT_32S,2)
#define MAT_32SC3 MAT_MAKETYPE(MAT_32S,3)
#define MAT_32SC4 MAT_MAKETYPE(MAT_32S,4)
#define MAT_32SC(n) MAT_MAKETYPE(MAT_32S,(n))

#define MAT_32FC1 MAT_MAKETYPE(MAT_32F,1)
#define MAT_32FC2 MAT_MAKETYPE(MAT_32F,2)
#define MAT_32FC3 MAT_MAKETYPE(MAT_32F,3)
#define MAT_32FC4 MAT_MAKETYPE(MAT_32F,4)
#define MAT_32FC(n) MAT_MAKETYPE(MAT_32F,(n))

#define MAT_64FC1 MAT_MAKETYPE(MAT_64F,1)
#define MAT_64FC2 MAT_MAKETYPE(MAT_64F,2)
#define MAT_64FC3 MAT_MAKETYPE(MAT_64F,3)
#define MAT_64FC4 MAT_MAKETYPE(MAT_64F,4)
#define MAT_64FC(n) MAT_MAKETYPE(MAT_64F,(n))

#define MAT_AUTO_STEP  0x7fffffff
//#define MAT_WHOLE_ARR  cvSlice( 0, 0x3fffffff )

#define CMAT_AUTO_STEP      0

#define MAT_MAT_CN_MASK          ((MAT_CN_MAX - 1) << MAT_CN_SHIFT)
#define MAT_MAT_CN(flags)        ((((flags) & MAT_MAT_CN_MASK) >> MAT_CN_SHIFT) + 1)
#define MAT_MAT_TYPE_MASK        (MAT_DEPTH_MAX*MAT_CN_MAX - 1)
#define MAT_MAT_TYPE(flags)      ((flags) & MAT_MAT_TYPE_MASK)
#define MAT_MAT_CONT_FLAG_SHIFT  14
#define MAT_MAT_CONT_FLAG        (1 << MAT_MAT_CONT_FLAG_SHIFT)
#define MAT_IS_MAT_CONT(flags)   ((flags) & MAT_MAT_CONT_FLAG)
#define MAT_IS_CONT_MAT          MAT_IS_MAT_CONT
//#define MAT_SUBMAT_FLAG_SHIFT    15
//#define MAT_SUBMAT_FLAG          (1 << MAT_SUBMAT_FLAG_SHIFT)
//#define MAT_IS_SUBMAT(flags)     ((flags) & MAT_MAT_SUBMAT_FLAG)

#define MAT_CONTINUOUS_FLAG		MAT_MAT_CONT_FLAG

#define MAT_ELEM_SIZE(type) \
	(MAT_MAT_CN(type) << ((((sizeof(size_t)/4+1)*16384|0x3a50) >> MAT_MAT_DEPTH(type)*2) & 3))

//#define MAT_MAGIC_MASK          0xFFFF0000
//#define MAT_MAT_MAGIC_VAL    0x42420000
//#define MAT_TYPE_NAME_MAT    "matrix"

typedef unsigned char uchar;

class CSize;
class CPoint;
class CRect;

//////////////////////////////////////////////////////////////////////////////
// CPoint is a 2D point, (x, y)
class CPoint
{
public:
	CPoint(); // default constructor, set x and y as 0
	CPoint(int _x, int _y); // x and y coordinate
	CPoint(const CPoint& pt); // copy constructor
	CPoint(const CSize& sz); // construct CPoint from CSize

	CPoint& operator = (const CPoint& pt);
	inline bool inside( const CRect& r ) const;

	int x, y; // x and y coordinate
};

class CPoint2f
{
public:
	CPoint2f(); // default constructor, set x and y as 0
	CPoint2f(float _x, float _y); // x and y coordinate
	CPoint2f(const CPoint2f& pt); // copy constructor
	//CPoint2f(const CSize& sz); // construct CPoint from CSize

	CPoint2f& operator = (const CPoint2f& pt);
	//inline bool inside(const CPoint2f& r) const;

	float x, y; // x and y coordinate
};

////////////////////////////////////////////////////////////////////////////////////
// CSize has w and h, used for matrix or image
class CSize
{
public:
	CSize(); // default constructor, set w and h as 0
	CSize(int _width, int _height);
	CSize(const CSize& sz); // copy constructor
	CSize(const CPoint& pt); // construct CSize from CPoint

	CSize& operator = (const CSize& sz);
	//! the area (width*height)
	int area() const;

	int width, height; // the width and the height
};

/////////////////////////////////////////////////////////////////////////////////////////
// CRect is a rectangle, with x, y coordinate of top-left corner, and w, h
class CRect
{
public:
	CRect();
	CRect(int _x, int _y, int _width, int _height);
	CRect(const CRect& r);
	CRect(const CPoint& org, const CSize& sz);
	CRect(const CPoint& pt1, const CPoint& pt2);

	CRect& operator = ( const CRect& r );
	//! the top-left corner
	CPoint tl() const;
	//! the bottom-right corner
	CPoint br() const;

	//! size (width, height) of the rectangle
	CSize size() const;
	//! area (width*height) of the rectangle
	int area() const;

	//! checks whether the rectangle contains the point
	bool contains(const CPoint& pt) const;

	int x, y, width, height; //< the top-left corner, as well as width and height of the rectangle
};

class CRect2f
{
public:
	CRect2f();
	CRect2f(float _x, float _y, float _width, float _height);
	CRect2f(const CRect2f& r);
	//CRect(const CPoint& org, const CSize& sz);
	//CRect(const CPoint& pt1, const CPoint& pt2);

	CRect2f& operator = (const CRect2f& r);
	////! the top-left corner
	//CPoint tl() const;
	////! the bottom-right corner
	//CPoint br() const;

	////! size (width, height) of the rectangle
	//CSize size() const;
	////! area (width*height) of the rectangle
	//int area() const;

	////! checks whether the rectangle contains the point
	//bool contains(const CPoint& pt) const;

	float x, y, width, height; //< the top-left corner, as well as width and height of the rectangle
};

//class  CRange
//{
//public:
//	CRange();
//	CRange(int _start, int _end);
//	int size() const;
//	bool empty() const;
//	static CRange all();
//
//	int start, end;
//};

///////////////////////////////////////////////////////////////////////////////////////
// CMat is a 2D matrix, which support multi-channel
// char, unsigned char, short, int, float, double are supported
// matrix has reference counter
class CMat
{
public:
    // default constructor
    CMat();

    // constructs 2D matrix of the specified size and type
    // (_type is MAT_8UC1, MAT_64FC3, MAT_32SC(12) etc.)
    CMat(int rows, int cols, int type);
    CMat(CSize size, int type);

    // copy constructor
    CMat(const CMat& m);

	CMat(int rows, int cols, int type, void* data, size_t step = CMAT_AUTO_STEP);
    CMat(CSize size, int type, void* data, size_t step = CMAT_AUTO_STEP);
    //CMat(int ndims, const int* sizes, int type, void* data, const size_t* steps = 0);
 
    // destructor - calls release()
    ~CMat();

    // assignment operators
    CMat& operator = (const CMat& m);

	//2017.4.16 jimmypan added
	CMat mul(CMat &m);
	//CMat operator * (float fval);
	//CMat operator + (float fval);	//-！如果复数矩阵，只加实部，不加虚部！
	//CMat operator + (CMat &m);

    // allocates new matrix data unless the matrix already has specified size and type.
    // previous data is unreferenced if needed.
    void create(int rows, int cols, int type);
    void create(CSize size, int type);
    void create(int ndims, const int* sizes, int type);

    // increases the reference counter; use with care to avoid memleaks
    void addref();

    // decreases reference counter;
    // deallocates the data when reference counter reaches 0.
    void release();

    // deallocates the matrix data
    void deallocate();

    // internal use function; properly re-allocates _size, _step arrays
    void copySize(const CMat& m);

	//size_t elemSize1() const;

    // returns element type, similar to MAT_MAT_TYPE(cvmat->type)
    int type() const;

    // returns element type, similar to MAT_MAT_DEPTH(cvmat->type)
    int depth() const;
    
	// returns element type, similar to MAT_MAT_CN(cvmat->type)
    int channels() const;
    
	// returns step/elemSize1()
    //size_t step1(int i=0) const;

    // returns true if matrix data is NULL
    bool empty() const;

    // returns the total number of matrix elements
    size_t total() const;

    // returns pointer to i0-th sub-matrix along the dimension #0
    uchar* ptr(int i0=0);
    const uchar* ptr(int i0=0) const;

    // returns pointer to (i0,i1) sub-matrix along the dimensions #0 and #1
    uchar* ptr(int i0, int i1);
    const uchar* ptr(int i0, int i1) const;

    // returns pointer to (i0,i1,i3) sub-matrix along the dimensions #0, #1, #2
    uchar* ptr(int i0, int i1, int i2);
    const uchar* ptr(int i0, int i1, int i2) const;

    // returns pointer to the matrix element
    uchar* ptr(const int* idx);

    // returns read-only pointer to the matrix element
    const uchar* ptr(const int* idx) const;

     //includes several bit-fields:
     //    - the magic signature
     //    - continuity flag
     //    - depth
     //    - number of channels
    int flags;

    // the matrix dimensionality, >= 2
    int dims;
    
	// the number of rows and columns or (-1, -1) when the matrix has more than 2 dimensions
    int rows, cols;

    // pointer to the matrix data
    uchar* data;

    // pointer to the reference counter;
    // when matrix points to user-allocated data, the pointer is NULL
    int* refcount;

    //! helper fields used in locateROI and adjustROI
    uchar* datastart;
    uchar* dataend;
    uchar* datalimit;
	
	// size of matrix for each dimension
    struct MSize
    {
        MSize(int* _p);
        CSize operator()() const;
        const int& operator[](int i) const;
        int& operator[](int i);
        operator const int*() const;
        bool operator == (const MSize& sz) const;
        bool operator != (const MSize& sz) const;

        int* p;
    };

	// step for each dimension
    struct MStep
    {
        MStep();
        MStep(size_t s);
        const size_t& operator[](int i) const;
        size_t& operator[](int i);
        operator size_t() const;
        MStep& operator = (size_t s);

        size_t* p;
        size_t buf[2];
    protected:
        MStep& operator = (const MStep&);
    };

    MSize size;
    MStep step;

protected:
    void initEmpty();
};
CMat operator * (CMat m, float fval);
CMat operator + (CMat m, float fval);
CMat operator + (CMat m, CMat n);

inline void printCMat2f(CMat& m)
{
	printf("+++++++++++++++ begin +++++++++++++++\n");
	int ch = m.channels();
	for (int i = 0; i < m.rows; i++)
	{
		float* pdata = (float*)m.ptr(i);
		for (int j = 0; j < m.cols; j++)
		{
			for (int k = 0; k < ch; k++)
			{
				printf("%f,", pdata[j*ch + k]);
			}
			printf(" ");
		}
		printf("\n");
	}
	printf("================= end =================\n");
}
inline void printCMat2i(CMat& m)
{
	printf("+================= begin =================\n");
	int ch = m.channels();
	for (int i = 0; i < m.rows; i++)
	{
		uchar* pdata = m.ptr(i);
		for (int j = 0; j < m.cols; j++)
		{
			for (int k = 0; k < ch; k++)
			{
				printf("%3d,", pdata[j*ch + k]);
			}
			printf(" ");
		}
		printf("\n");
	}
	printf("================= end =================\n");
}
#endif // end of the file
