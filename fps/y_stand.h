/********************************************************************/
/* y_stand.h is the declaration part for fast program               */
/********************************************************************/
#ifndef Y_STAND_HDR						/* Make sure that the header file is expanded in */
#define Y_STAND_HDR						/* every .c file at most once.					 */


#include <windows.h>


	#include <stdio.h>
	#include <math.h>
	#include <stdlib.h>


#define and 	&&			/* some local defines */
#define or  	||
#define mod 	%
#define not 	!
#define is_not 	!=
#define is	==
#define PI	3.141592654

//#define InSize(N)     (((N)>=0) and ((N)<PicSize))
//#define InFrame(I,J)  (InSize(I) and InSize(J))
#define Sqr(x)  ((x)*(x))
#define Sqrt(x) (sqrt(x))

//typedef enum boolean{FALSE = 0, TRUE = 1} boolean;

/*big real positive number */
#define MAX_DISTANCE  	9999999    

#define MaxTri	 32	/* The maximum number of triangles that a vertex can participate in. */

#define ABS(x) 		((x)>0?(x):(-(x)))
#define MAX(x,y) 	((x)>(y)?(x):(y))
#define MIN(x,y) 	((x)>(y)?(y):(x))
#define FAST_MIN(x,y) (((x) > (y)) ? (x) = (y) : 0)	/* Can replace a macro call of the form
													   x = MIN(x, y)					*/
#define MIN4(a,b,c,d) 	MIN(MIN(a,b),MIN(c,d))
//#define MAX4(a,b,c,d) 	MAX(MAX(a,b),MAX(c,d))

#define Alive 0		/* Alive and Far are the values that are placed in the Back Pointer array when the vertex is not */
#define Far   -1	/* on the heap.																					 */


// allocate 2-D matrix: [hsize][vsize] of type _MTYPE
#define matrix_allocate(matrix,_vsize,_hsize,_MTYPE) { \
	int _i; \
	_MTYPE *myTypePtr=NULL; \
	matrix=(_MTYPE **)malloc(_vsize*sizeof(myTypePtr)); \
	if (matrix==NULL) { \
		fprintf(stderr,"Error: Can't allocate memory for 2D array\n"); \
		exit(1); \
	} \
	for (_i=0 ; _i<_vsize; _i++) { \
		matrix[_i]=(_MTYPE *)malloc(_hsize*sizeof(_MTYPE)); \
		if (matrix[_i] ==NULL) { \
			fprintf(stderr,"Error: Can't allocate memory for 2D array\n"); \
			exit(1); \
		} \
	} \
}

// free 2-D matrix: with _vsize lines
#define matrix_free(matrix,_vsize) { \
	int _i; \
	for (_i=0 ; _i<_vsize; _i++) \
		free(matrix[_i]); \
	free(matrix); \
}



struct heap_element { 		/* element of the heap													*/
	double		u;			/* element value (u)													*/
	int 		v; 			/* back pointer to the vertex (the index of the vertex in the array V)	*/
};


class Heap {

	public:

		Heap(int Vnum) { init_heap(Vnum); }
		~Heap() { free_heap(); }

		void reset();
		void upheap(int k);
		void insert(double u, int v);
		void downheap(int k);
		void remove_top();

	protected:

		void init_heap(int Vnum);
		void free_heap(void);

	public:

		int			*BP;			/* back pointer to place in the list */
		int			N;				/* number of elements in the heap array */
		heap_element *a; 			/* heap array */

};




struct Triangle { 	/* defines one triangle					*/
	int	Vind[3];	/* Index for the 3 vertices				*/
	double  b[5];   /* surface cooefficients				*/
					/* Du = (2b0x+b2y+b3,2b1y+b2x+b4)		*/
					/* for Vind[0]-Vind[1] = x exis			*/
	boolean Visited;/* Determine whether this vertex was
					   already visited when backtracking to
					   find a smooth geodesic.				*/
	boolean Split;	/* This field insdicates whether on of
					   the edges of this triangle is
					   splitting edge.						*/
};

struct Stencil {		/* virtual numerical connections  */
				/*   true triangles in most cases */
	double	Ctheta;		/* cos(theta) angle between the edges */
	double  Stheta; /* sin(theta)^2 when theta is the angle between the edges */
	int	v1,v2;		/* index of neigboring vertices */
	double	l1,l2;		/* edge length to v1 and to v2 */
	double Ctheta_mul_l1_div_l2_minus_one_mul_2, Ctheta_mul_l2_div_l1_minus_one_mul_2;
								/* Remember the value of the l1 and l2 multiplies
										    by Ctheta minus one	multiplied by 2	*/
	double sqr_l1, sqr_l2; /* The values of l1 and l2 multiplied by themselves. */
						/* Replaces 1 - a2 * (ctheta_mul_mb_div_ma_minus_one_mul_2 + 1) / b2 */ 
	double shortcut1_1, shortcut1_2;
						/* Replaces - Stheta * a2*/ 
	double shortcut2_1, shortcut2_2;

};

struct Vertex { 		/* defines one Numerical Graph vertex */
	double	x,y,z;		/*x,y,z coordinates   */
	double	U;	   	/* U vlaue	*/
	int	Tind[MaxTri];	/* link back to triangles */
					/* MaxTri = Max # triangles at one vertex */
	int	si,vn;		/* number of vertex connections */
					/* si is the #of ST, vn is the index to VN*/
	int ti;			/* The number of triangles this vertex is member of */
	struct  Stencil	ST[MaxTri];	/* numerical connection, 
						updating stenciles  */
	int	VN[3*MaxTri]; /* Neighboring dirrectional vertices indexes*/
				/* Vertex to update */
	int	ST_For_VN[3*MaxTri]; /* The first stenciel of the neighbour with this vertex in it */

	int Split[2];	/* The indexes of the vertexes that split the triangles of that vertex
					   or NO_DATA if there is no vertex.									*/
							 
	boolean source;	  /* Indicates wether the variable is a source in one of the marches
					     (meaning it appears in the xStart, yStart variables).				*/
	boolean current_source;	/* Is it a source that got the U value 0 in the last run of
						 the fast march algoritem.											*/
	double MaxU;
};

#define NoData -1					/* Indicate no data for this field.						*/


struct Point {						/* defines a point in 3D with a U value					*/
	double  x,y,z;					/* x,y,z coordinates IN 3D								*/
	double  U;						/* The U value of the point.							*/
};

struct Vector {						/* defines a vector in 3D								*/
		double  x,y,z;				/* x,y,z coordinates of the vector in 3D.				*/
};


#endif /* end of #ifndef Y_STAND_HDR */
