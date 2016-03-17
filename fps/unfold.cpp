/***************************************************************************/
/* Triangle.c : test fast marching Eikonal on acute (non obtuse)           */ 
/* triangulated grids	   											  	   */
/* This program was written by Ronny Kimmel    						       */
/***************************************************************************/


#include "y_stand.h"									/* declaration file */
#include "mex.h"
#include <stdlib.h>
#include "string.h" 


class FMM {

public:

	FMM(const double *X, const double *Y, const double *Z, const int *TRI, int vnum, int tnum) : Vnum(vnum), Tnum(3*tnum) {
		CombineVectors(&V, &T, &NonSplitTnum, 
					X,Y,Z,TRI, Vnum, Tnum); 
        InitGraph(T,V, Vnum, NonSplitTnum, Tnum); /* init V to isolated vertices */
		heap = new Heap(Vnum);
    }
	void March(double *SourceVal, double *DistanceData) {
        heap->reset();
        MarchAlg(V, Vnum, SourceVal); //, SourceVal+Vnum);	/* Perform the march */			
		for(int i = 0; i < Vnum; i++){	/* Fill up the result matrix with data.					*/
           *(DistanceData++) = V[i].U;
        }
    }

	~FMM() {
        /* Free the resources used by the program. */
		delete heap;
        free(T);
        free(V);
    }

	int GetV() { return Vnum; }
	int GetT() { return Tnum; }

protected:

	void InitGraph(struct Triangle * T, struct Vertex * V, int Vnum, int NonSplitTnum, int Tnum);
	void CombineVectors(struct Vertex *(V[]), struct Triangle *(T[]), int *NonSplitTnum, 
						const double *X, const double *Y, const double *Z, const int *TrianglesData, 
						int Vnum, int Tnum);
	void InitMarch(double srcval[], int Vnum);
	void MarchAlg(struct Vertex  *temp_V,	
				  int temp_Vnum,
				  double *srcval);
	void __forceinline __fastcall FMM::CallTUpdateNeighbors(int             i,	/* vertex number */
															int start_from		/* The source from which the current march began.	*/
															);
	__forceinline void TUpdateNeighbors(
		int    i,						/* The vertex to update				   */
		int	   becomes_alive			/* The vertex that now becomes alive.  */
	);
	__forceinline double __fastcall update(struct Stencil *current_stencil, int k, int l);

private:

	struct Triangle *T; 	/* triangles, numbered 1..Tnum */
	struct Vertex   *V; 	/*  vertices coordinates 1..VNum */
	int	Tnum, Vnum;
	int NonSplitTnum;		/* The number of triangles with no splitting. */
							/* Names of files that are used for writing the results. */
	Heap *heap;

};

/***************************************************************************/
/* L2  error the diff between the distance and U*/
/***************************************************************************/
double L2(Vertex *V, int Vnum)
{
	int             i;
	double          d, x, y, z;
	double          sum = 0, dxdy = 1.0/((double)Vnum);
	/* find the zero point */
	for (i = 0; i < Vnum; i++)
		if (V[i].U < 0.00000001) {
			x = V[i].x;
			y = V[i].y;
			z = V[i].z;
			mexPrintf ("Source at: (%g,%g,%g)\n",x,y,z);
		}
	for (i = 0; i < Vnum; i++) {
		d = sqrt(Sqr(x - V[i].x) + Sqr(y - V[i].y) + Sqr(z - V[i].z));
		sum += Sqr(V[i].U - d)*dxdy; 
	}
	return (sqrt(sum));
}
/***************************************************************************/
/* L1  error norm the diff between the distance and U*/
/***************************************************************************/
double L1(Vertex *V, int Vnum)
{
	int             i;
	double          d, x, y, z;
	double          sum = 0.0, dxdy = 1.0/((double)Vnum);
	/* find the zero point */
	for (i = 0; i < Vnum; i++)
		if (V[i].U < 0.00000001) {
			x = V[i].x;
			y = V[i].y;
			z = V[i].z;
		}
	for (i = 0; i < Vnum; i++) {
		d = sqrt(Sqr(x - V[i].x) + Sqr(y - V[i].y) + Sqr(z - V[i].z));
		sum += ABS(V[i].U - d)*dxdy;
	}
	return (sum);
}

/***************************************************************************/
/* CosAngle the cos of the angle at the vertex v0, between v1 and v2	   */
/***************************************************************************/
double CosAngle(int v0,int v1,int v2,struct Vertex *V)
{
	double x1,x2,y1,y2,z1,z2,res;
	if(v0 != -1 and v1 != -1 and v2 != -1){
		x1 = V[v1].x - V[v0].x;
		x2 = V[v2].x - V[v0].x;
		y1 = V[v1].y - V[v0].y;
		y2 = V[v2].y - V[v0].y;
		z1 = V[v1].z - V[v0].z;
		z2 = V[v2].z - V[v0].z;
		res = x1*x2+y1*y2+z1*z2;		/* dot product */
		res /= sqrt(x1*x1+y1*y1+z1*z1); /* normalize */
		res /= sqrt(x2*x2+y2*y2+z2*z2);
		return(res);
	}
	else
		return 0;
}
/***************************************************************************/
/* Length between the vertex v0 and v1									   */
/***************************************************************************/
double
Length(int v0,int v1,struct Vertex *V)
{
	double x1,y1,z1,res;
	if(v0 != -1 and v1 != -1){
		x1 = V[v1].x - V[v0].x;
		y1 = V[v1].y - V[v0].y;
		z1 = V[v1].z - V[v0].z;
		res = sqrt(x1*x1+y1*y1+z1*z1); /* distance */
		return(res);
	}
	else
		return MAX_DISTANCE;
}
/***************************************************************************/
/* nextT next triangle to be unfolded. find the triangle #, and the vertex */
/* Returns true is the next triangle was found.							   */
/* v1 and v2 indicate the edge that is common to the original triangle	   */
/* and the triangle to be unfolded.										   */
/* ti is the original triangle and v3 is the other vertex of the triangle  */
/* to be unfolded.														   */
/* vn is the index of the triangle to be unfolded.						   */
/***************************************************************************/
boolean
nextT(int ti, int v1, int v2, int *v3, int *tn,
      struct Triangle * T, struct Vertex * V, int Tnum)
{
	boolean				found = FALSE;	/* Indicates whether we found the next triangle. */
	int					i,				/* Index for the loop.							 */
						tj,				/* A candidate tp be the next triangle.			 */
						k;				/* Index for the inner loop.					 */
	/* scan every triangle of vi */
	for (i = 0; i < V[v1].ti and not found; i++) {
		tj = V[v1].Tind[i];
		if (tj < Tnum and tj != ti)
			/* search for tj in the list of v2 */
			for (k = 0; k < V[v2].ti and not found; k++) 
				if (V[v2].Tind[k] == tj && !T[tj].Split) {
					found = TRUE;
					*tn = tj;
				}
	}
	if (found){ /* find v3, the other vertex in the triangle to be unfolded.			 */
		if(T[*tn].Vind[0] == v1){
			if(T[*tn].Vind[1] == v2)
				*v3 = T[*tn].Vind[2];
			else
				*v3 = T[*tn].Vind[1];
		}
		else if(T[*tn].Vind[1] == v1){
			if(T[*tn].Vind[0] == v2)
				*v3 = T[*tn].Vind[2];
			else
				*v3 = T[*tn].Vind[0];
		}
		else{
			if(T[*tn].Vind[0] == v2)
				*v3 = T[*tn].Vind[1];
			else
				*v3 = T[*tn].Vind[0];
		}
	}
	return (found);
}

/***************************************************************************/
/* Split obtuse angles by unfolding splitting and connecting, return the   */
/* number of unfoldings that were nessesery.							   */
/* ti is the tirnalge to be splittined, V0 is the vertex with the obtuse   */
/* angle while V1 and V2 are the other vertexes of ti.					   */
/***************************************************************************/
int
Split(int ti, int V0, int V1, int V2, struct Triangle * T, struct Vertex * V,
      int NonSplitTnum, int Vnum)
{
	double          xv1,x1, y1, yv1,
					x2,				/* The distance between V0 and V2 */
					y2,
					x3,y3, xt2, xt3, yt3,
					e0,				/* The distance between v1 and v2 */
					e1,				/* The distance between v2 and v3.*/
					e2,				/* The distance between v0 and v1 */
					ta,				/* Tan of alpha.				  */
					cb, sb;			/* Cos and Sin of beta.			  */
	int             v1 = V1, v2 = V2, v3,
					tm=ti,			/* The current triangle we are
									   working on.					  */
					tn,				/* The triangle returned by NextT */
					count = 0;		/* The number of triangles unfolded
									   so far.						  */
									/* Becomes true when the split was done */
	boolean         splitter = FALSE;
	x2 = Length(V0, V2, V);
	y2 = 0;
	e0 = Length(V1, V2, V);
	e2 =  Length(V0, V1, V);
	xv1 = x1 = (x2 * x2 + e2 * e2 - e0 * e0) / (2.0 * x2);/* translation */
	yv1 = y1 = sqrt(e2 * e2 - x1 * x1);
	ta = -x1 / y1;		/* tan (alpha) in Fig. 1 */
	/* if there is a next triangle and not splited */
	while (nextT(tm, v1, v2, &v3, &tn, T, V, NonSplitTnum) and (not splitter) ) {
		count++;
		tm = tn;		/* Update the wording triangle. */
		cb = (x2 - x1) / sqrt(Sqr(x2 - x1) + Sqr(y2 - y1));	/* cos beta */
		sb = sqrt(1 - cb * cb);								/* sin beta */
		if (y2 < y1)	/* Adjast the sign of SIN(beta).				*/
			sb *= -1.0;
		xt2 = Length(v1, v2, V);
		e1 = Length(v2, v3, V);
		e2 = Length(v1, v3, V);
		xt3 = (xt2 * xt2 + e2 * e2 - e1 * e1) / (2.0 * xt2);
		yt3 = sqrt(e2 * e2 - xt3 * xt3);
		x3 = cb * xt3 - sb * yt3 + x1;
		y3 = sb * xt3 + cb * yt3 + y1;

		if (x3 > 0 and y3/x3 > ta) {		/* if we found a splitter */
			splitter = TRUE;
											/* Add the stencils involving the
											   splitting edge.			*/
			V[V0].ST[V[V0].si].Ctheta = (x3*xv1+y3*yv1)/
				sqrt((xv1*xv1+yv1*yv1)*(x3*x3+y3*y3));
			V[V0].ST[V[V0].si].v1 = V1;
			V[V0].ST[V[V0].si].v2 = v3;
			V[V0].ST[V[V0].si].l1 = Length(V0, V1, V);
			if(V[V0].si == MaxTri - 1) {
			//	mexPrintf("Warning, too many triangles to one vertex, result quality will suffer\n");
			}
			else
				V[V0].ST[V[V0].si++].l2 = sqrt(x3*x3+y3*y3);


			V[V0].ST[V[V0].si].Ctheta = x3/sqrt(x3*x3+y3*y3);
			V[V0].ST[V[V0].si].v1 = v3;
			V[V0].ST[V[V0].si].v2 = V2;
			V[V0].ST[V[V0].si].l1 =sqrt(x3*x3+y3*y3);
			if(V[V0].si == MaxTri - 1) {
			//	mexPrintf("Warning, too many triangles to one vertex, result quality will suffer\n");
			}
			else
				V[V0].ST[V[V0].si++].l2 =Length(V0, V2, V);

			if(V[v3].vn == 3 * MaxTri - 1) {
			//	mexPrintf("Warning, too many triangles to one vertex, result quality will suffer\n");
			}
			else
				V[v3].VN[V[v3].vn++] = V0; /* add dirrectional edge
											  to v3					*/

			T[NonSplitTnum + (ti * 2)].Vind[0] = V1;	/* Add the triangles of the splitting. */
			T[NonSplitTnum + (ti * 2)].Vind[1] = V0;
			T[NonSplitTnum + (ti * 2)].Vind[2] = v3;
			T[NonSplitTnum + (ti * 2)].Split = TRUE;
			if(V[V1].ti == MaxTri - 1) {
			//	mexPrintf("Warning, too many triangles to one vertex, result quality will suffer\n");
			}
			else
				V[V1].Tind[V[V1].ti++] = NonSplitTnum + (ti * 2);
			if(V[V0].ti == MaxTri - 1) {
			//	mexPrintf("Warning, too many triangles to one vertex, result quality will suffer\n");
			}
			else
				V[V0].Tind[V[V0].ti++] = NonSplitTnum + (ti * 2);
			if(V[v3].ti == MaxTri - 1){
			//	mexPrintf("Warning, too many triangles to one vertex, result quality will suffer\n");
			}
			else
				V[v3].Tind[V[v3].ti++] = NonSplitTnum + (ti * 2);
			T[NonSplitTnum + (ti * 2) + 1].Vind[0] = v3;
			T[NonSplitTnum + (ti * 2) + 1].Vind[1] = V0;
			T[NonSplitTnum + (ti * 2) + 1].Vind[2] = V2;
			T[NonSplitTnum + (ti * 2) + 1].Split = TRUE;
			if(V[v3].ti == MaxTri - 1) {
			//	mexPrintf("Warning, too many triangles to one vertex, result quality will suffer\n");
			}
			else
				V[v3].Tind[V[v3].ti++] = NonSplitTnum + (ti * 2) + 1;
			if(V[V0].ti == MaxTri - 1) {
			//	mexPrintf("Warning, too many triangles to one vertex, result quality will suffer\n");
			}
			else
				V[V0].Tind[V[V0].ti++] = NonSplitTnum + (ti * 2) + 1;
			if(V[V2].ti == MaxTri - 1) {
			//	mexPrintf("Warning, too many triangles to one vertex, result quality will suffer\n");
			}
			else
				V[V2].Tind[V[V2].ti++] = NonSplitTnum + (ti * 2) + 1;


		}
		else {						   /* we have not found a splitter,
										  continue unfolding			*/
			if (x3 < 0){
				v1 = v3; x1 = x3; y1 = y3;
			} else {
				v2 = v3; x2 = x3; y2 = y3;
			}
		}	
	}
	return(count);				/* Return the number of triangles that were
								   unfolded.							   */
}

/***************************************************************************/
/* InitGraph init the numerical graph V, unfold in order to split obtuse   */
/* triangles.															   */
/***************************************************************************/
void
FMM::InitGraph(struct Triangle * T, struct Vertex * V, int Vnum, int NonSplitTnum, int Tnum)
{
	int             i, k,			/* Indexes for the loops.				*/
					ti,				/* Index for the current triangle.		*/
					v1,v2,			/* Indexes for the neighbours in the
									   triangle of the the current vertex.	*/
					count= 0,		/* The number of unfolding in one split.*/
					mcount=0;		/* The maximum value count got.			*/
	boolean         found;			/* Used for adding the neighbours to
									   every triangle.						*/
	double          ca;				/* The cosin value that determine
									   whether this triangle is obtuse.		*/
	struct Stencil *p_st;			/* Pointer to a stencil.				*/
	int ind, si_count;				/* Used for precalculting values for the
									   stencil.								*/


	/* Initialize the vertixes. */
	for (i = 0; i < Vnum; i++) {	/* zero counters of all vertices */
		V[i].si = 0;				/* # of connections to other triangles */
		V[i].vn = 0;				/* # of connections to other vertices */
	}

	/* Set the split field of the triangles that exist now, before the splitting to false.	*/
	for(i = 0; i < NonSplitTnum; i++)
		T[i].Split = FALSE;

	
	for (i = 0; i < Vnum; i++) {		/* scan all vertices */
		for (ti = 0; ti < V[i].ti; ti++){/* scan connected triangles */
			if (V[i].Tind[ti] < NonSplitTnum) {	/* if valid triangle */
												/* Make v1 and v2 the neighbours.			*/
				if (T[V[i].Tind[ti]].Vind[0] == i){
					v1 = T[V[i].Tind[ti]].Vind[1];
					v2 = T[V[i].Tind[ti]].Vind[2];
				}
				else if (T[V[i].Tind[ti]].Vind[1] == i){
					v1 = T[V[i].Tind[ti]].Vind[2];
					v2 = T[V[i].Tind[ti]].Vind[0];
				}
				else if (T[V[i].Tind[ti]].Vind[2] == i){
					v1 = T[V[i].Tind[ti]].Vind[0];
					v2 = T[V[i].Tind[ti]].Vind[1];
				}

				found = FALSE;					/* Add v1 as a neighbour if it is not already
												   a neighbour.								*/
				for (k = 0; k < V[i].vn; k++)
					if (v1 == V[i].VN[k])
						found = TRUE;
				if (not found)
					V[i].VN[V[i].vn++] = v1;

				found = FALSE;					/* Add v2 as a neigbour if it is not already
												   a neighbour.								*/
				for (k = 0; k < V[i].vn; k++)
					if (v2 == V[i].VN[k])
						found = TRUE;
				if (not found)
					V[i].VN[V[i].vn++] = v2;
			
				ca = CosAngle(i,v1,v2,V);
				if (ca < 0){					/* If this triangle is an obtuse angle		*/
					count = Split(V[i].Tind[ti],i,v1,v2,
						T,V,NonSplitTnum,Vnum);
					if (count > mcount)			/* Update m count.							*/
						mcount = count;
				} 
				else {							/* If no splitting was nessesery create
												   the stencil for this vertex and triangle.*/
					V[i].ST[V[i].si].Ctheta = ca;
					V[i].ST[V[i].si].v1 = v1;
					V[i].ST[V[i].si].l1 = Length(i,v1,V);
					V[i].ST[V[i].si].v2 = v2;
					V[i].ST[V[i].si].l2 = Length(i,v2,V);
					V[i].si++;
				}
			}
		}
	}

	for(ind = 0; ind < Vnum; ind++)					/* Calculate the data for each stencil.	*/
		for(p_st = V[ind].ST, si_count = V[ind].si - 1; si_count >= 0; p_st++, si_count--){
			p_st->Stheta =				1 - Sqr(p_st->Ctheta);
			p_st->Ctheta_mul_l1_div_l2_minus_one_mul_2 =p_st->Ctheta * p_st->l1 * 2 / p_st->l2 - 2;
			p_st->Ctheta_mul_l2_div_l1_minus_one_mul_2 =p_st->Ctheta * p_st->l2 * 2 / p_st->l1 - 2;
			p_st->sqr_l1 = p_st->l1 * p_st->l1;
			p_st->sqr_l2 = p_st->l2 * p_st->l2;
			p_st->shortcut1_1 = 1 - p_st->sqr_l1 * (p_st->Ctheta_mul_l2_div_l1_minus_one_mul_2 + 1) / p_st->sqr_l2;
			p_st->shortcut1_2 = 1 - p_st->sqr_l2 * (p_st->Ctheta_mul_l1_div_l2_minus_one_mul_2 + 1) / p_st->sqr_l1;
			p_st->shortcut2_1 = - p_st->Stheta * p_st->sqr_l1;
			p_st->shortcut2_2 = - p_st->Stheta * p_st->sqr_l2;
		}

//	mexPrintf("\nNumber of unfoldings = %d\n",mcount);
}


void FMM::CombineVectors(struct Vertex *(V[]), struct Triangle *(T[]), int *NonSplitTnum, 
					const double *X, const double *Y, const double *Z, const int *TrianglesData, 
					int Vnum, int Tnum){

	int i, j;										/* Indexes for the loops.				*/
	const double *CoordinatesData[3] = {X,Y,Z};		/* Pointer to the data in the three 
													   vector arrays.						*/
	
	*NonSplitTnum = Tnum / 3;

													/* Allocate memory for both triangles and
												       vertixes.							*/
    *T = (struct Triangle *) malloc(sizeof(struct Triangle) * Tnum);
    if(*T == NULL){
		fprintf(stderr, "Out of memory for triangles - exiting.\n");
		exit(-1);
	}
	*V = (struct Vertex *)   malloc(sizeof(struct Vertex) * Vnum);
    if(*V == NULL){
		free(T);
		fprintf(stderr, "Out of memory for vertiexes - exiting.\n");
		exit(-1);
	}
	
	for(i = 0; i < Vnum; i++){						/* Move the data to V and T.			*/
		(*V)[i].x = *((double *)CoordinatesData[0]++);
		(*V)[i].y = *((double *)CoordinatesData[1]++);
		(*V)[i].z = *((double *)CoordinatesData[2]++);
	}
	for(i = 0; i < 3; i++)
		for(j = 0; j < *NonSplitTnum; j++){
			(*T)[j].Vind[i] = *((int *)TrianglesData++);
		}
		
	for(i = 0; i < Vnum; i++){						/* Add every triangle to its vertixes.	*/
		(*V)[i].ti = 0;
	}
	/* Can be greatly improved! */
	for(i = 0; i < *NonSplitTnum; i++){
		if((*T)[i].Vind[0] != *NonSplitTnum){
			(*V)[(*T)[i].Vind[0]].Tind[(*V)[(*T)[i].Vind[0]].ti++] = i;
			(*V)[(*T)[i].Vind[1]].Tind[(*V)[(*T)[i].Vind[1]].ti++] = i;
			(*V)[(*T)[i].Vind[2]].Tind[(*V)[(*T)[i].Vind[2]].ti++] = i;
		}
	}
	
	return;

}




/***************************************************************************/
/* IntiMarch procedure, Init phase of the procedure, the march is from the */
/* points specified by xplace and yplace and the index start_from that	   */
/* the algoritem gets. If the flag ALL_SOURCES_TOGETHER is defined the	   */
/* the procedure ignore the third argument and use all points in xplace	   */
/* and yplace as sources.												   */
/***************************************************************************/
void FMM::InitMarch(double srcval[], int Vnum){
	int i;
	for(i = 0; i < Vnum; i++){

		if (srcval[i] >= MAX_DISTANCE) {
			V[i].U = MAX_DISTANCE;					/* Initialize the distance				*/
			heap->BP[i] = Far;							/* Originally all are far.				*/
			V[i].source = FALSE;
			//V[i].source_num = 0;
		}
		else {
			V[i].U = srcval[i];
			heap->insert(V[i].U, i);
			V[i].source = TRUE;
			//V[i].source_num = srcnum[i];	/* Set the source num field.			*/
		}

	}


	return;

}

/***************************************************************************/
/* Quadratic.c  solving a quadratic equation qa x^2+qb x+qc=0              */
/* Only the + solution is returned.										   */
/***************************************************************************/
__forceinline boolean __fastcall
Quadratic(double qa, double qb, double qc, double *x1)
{

	double d;

	d = qb*qb - (qa+qa)*(qc+qc); /* Discremenat */

    if ( d >= 0 ){  /* in case the Discremenat >= 0 */
		*x1 = (sqrt(d) - qb)/(qa + qa);
        return TRUE;
	}
    else {
        return FALSE;
    }
}

/***************************************************************************/
/* update find the solution to one triangle problem                        */
/***************************************************************************/
__forceinline double __fastcall
FMM::update(struct Stencil *current_stencil, int k, int l){

	double u1, u2;
	double u;
	double t;
	double /*ctheta_mul_mb_div_ma_minus_one_mul_2,*/ ctheta_mul_ma_div_mb_minus_one_mul_2;
																			/* Both of them are between 0 and 1 becuase
																			   the triangles are not obuse.				*/

	u1 = V[k].U;
	u2 = V[l].U;
	u = u2 - u1;
	if (u >= 0) {
		ctheta_mul_ma_div_mb_minus_one_mul_2 = current_stencil->Ctheta_mul_l2_div_l1_minus_one_mul_2;
																			/* A shortcut */
		if (Quadratic(current_stencil->shortcut1_2,							/* If the quadratic equation has solutions */
			      u * ctheta_mul_ma_div_mb_minus_one_mul_2,
			      (u * u + current_stencil->shortcut2_2), &t)
				  and (-(u + u) >= ctheta_mul_ma_div_mb_minus_one_mul_2 * t)
				  and current_stencil->Ctheta_mul_l1_div_l2_minus_one_mul_2*(t-u) <= (u + u))
				return (u1 + t);
		else									/* If the quadratic equation has no solution,
												   or the solutions are out of the trinagle */
			return MAX_DISTANCE;
	}
	else{
		ctheta_mul_ma_div_mb_minus_one_mul_2 = current_stencil->Ctheta_mul_l1_div_l2_minus_one_mul_2;
																		/* A shortcut */

		if (Quadratic(current_stencil->shortcut1_1,						/* If the quadratic equation has solutions */
			      -u * ctheta_mul_ma_div_mb_minus_one_mul_2,
			      (u * u + current_stencil->shortcut2_1), &t)
				  and (u + u >= ctheta_mul_ma_div_mb_minus_one_mul_2 * t)
				  and current_stencil->Ctheta_mul_l2_div_l1_minus_one_mul_2*(t+u) <= -(u + u))
				return (u2 + t);
		else									/* If the quadratic equation has no solution,
												   or the solutions are out of the trinagle */
			return MAX_DISTANCE;
	}


	
}

/***************************************************************************/
/* Update Neighbors procedure, update the U values of all the neighbours of*/
/* a vertex in a triangulation.											   */ 
/***************************************************************************/
__forceinline void
FMM::TUpdateNeighbors(
	int    i,						/* The vertex to update				   */
	int	   becomes_alive			/* The vertex that now becomes alive.  */
){
	double          u, t3, u0;
	int             n = heap->BP[i];
	int             k, l;
	int				st_ind;			/* Index for the stencils.				*/
	struct Stencil	*st_ptr;		/* Pointer to the current stencil.		*/


	u = V[i].U;
	u0 = u;
														/* Do for every stencil.						*/
	for (st_ptr = V[i].ST, st_ind = V[i].si; st_ind > 0; st_ind--, st_ptr++) {
														/* check all numerical connections (triangles)	*/
		k = st_ptr->v1;
		l = st_ptr->v2;
		if(k == becomes_alive || l == becomes_alive){	/* Use this stencil to update only if one k or
														   l is the new alive vertex.					*/
			if(k == becomes_alive)						/* Do the one dimensional update.				*/
				u = MIN(u, V[k].U + st_ptr->l1);
			else
				u = MIN(u, V[l].U + st_ptr->l2);
			if(heap->BP[k] == Alive || heap->BP[l] == Alive){		/* Do update only if the two other vertexes
														   of the stencil are alive, otherwise this
														   stencil will be called again anyway.		*/
				t3 = update(st_ptr, k, l);

				FAST_MIN(u, t3);						/* Minimize the distance.			*/
			}
			

//			if (u0 > u){
//				V[i].source_num = MAX(V[l].source_num, V[k].source_num);
//			}

		}
	}

//	if (V[i].source) { u = V[i].U; }

	if (n == Far){		/* not in heap                 */
		heap->insert(u, i);	/* also changes BP to Trail	   */
		V[i].U = u;
	}
	else {				/* change position in the heap */
		if (heap->a[n].u > u) {
			heap->a[n].u = V[i].U = u;
			heap->upheap(n);
		}
	}
}
/***************************************************************************/
/* Update Triangle Neighbors procedure					   */
/***************************************************************************/
void __forceinline
__fastcall FMM::CallTUpdateNeighbors(
	int             i,	/* vertex number */
	int start_from		/* The source from which the current march began.	*/
){
	int local_vn_count;			/* The number of neighbours yet to cover. */
	int *VN_ind;					/* Pointer to the array of neighbours */
										/* Call TUpdateNeighbors for every neighbour vertex of
										   the given vertex. */

//	mexPrintf("%d  (%d)\n", i,start_from);

	for (VN_ind = V[i].VN, local_vn_count = 0; local_vn_count < V[i].vn; VN_ind++, local_vn_count++)
										/* If the neighbour is not alive					*/
		if(heap->BP[*VN_ind] is_not Alive /*&& V[*VN_ind].source == FALSE*/){
			TUpdateNeighbors(*VN_ind, i);
		}

}

/***************************************************************************/
/* MarchAlg procedure, for fast marching method, extentions of Sethian     */
/* This procedure works with triangulations that contains only non-obtuse  */
/* triangles															   */
/***************************************************************************/
void FMM::MarchAlg(
	struct Vertex  *temp_V,	/* vertices coordinates 1..Vnum */
	int             temp_Vnum,
	double *srcval){
	V = temp_V;						/* Save the array of vertexes in the V array.		*/
	Vnum = temp_Vnum;				/* Save the number of vertexes in the array V.		*/

	InitMarch(srcval, Vnum);

									/* CAllTUpdateNeighbours for every vertex when it becomes
									   the smallest vertex in the heap. Exit the loop when
									   there are no more sources which are not alive */
	


	while (heap->N != 0){
		CallTUpdateNeighbors(heap->a[1].v, 0);
		heap->remove_top();
	}


}






/***************************************************************************/
/* This function is the entry point of the dll, matlab shuld call this	   */
/* function with a surface matrix and a sources matrix, the dll return the */
/* distances matrix and can print various types of output, according to	   */
/* the compilation flags that are defined.								   */
/***************************************************************************/
void mexFunction(
    int           nlhs,           /* number of expected outputs */
    mxArray       *plhs[],        /* array of pointers to output arguments */
    int           nrhs,           /* number of inputs */
    const mxArray *prhs[]         /* array of pointers to input arguments */
)
 {

    int i,j;
	int dims[2];


    char command[1024];
	int strlen = mxGetN(prhs[0]);
	strlen = strlen < mxGetM(prhs[0]) ? mxGetN(prhs[0]) : strlen;
	strlen += 1;
	mxGetString (prhs[0], command, strlen);
    
    prhs++;
    nrhs--;

  	if (!strcmpi(command, "init")) {
		FMM *fmm = new FMM( (double *)mxGetData(prhs[1]), (double *)mxGetData(prhs[2]), (double *)mxGetData(prhs[3]), (int *)mxGetData(prhs[0]),
			           MAX(mxGetM(prhs[1]), mxGetN(prhs[1])), mxGetM(prhs[0]) );        

        dims[0] = 1;
        dims[1] = 1;
        plhs[0] = mxCreateNumericArray(2, dims, mxUINT32_CLASS, mxREAL);
        ((unsigned int *)mxGetPr(plhs[0]))[0] = (unsigned int)fmm;

    }
    else if (!strcmpi(command, "march")) {

	    FMM *fmm = (FMM*)((unsigned int *)mxGetPr(prhs[0]))[0];

    	int Snum = mxGetN(prhs[1]); // Number of sources
		int Vnum = fmm->GetV();

        plhs[0] = mxCreateDoubleMatrix(Vnum, Snum, mxREAL);

	    double *SourceVal;

        if (mxGetM(prhs[1]) == 1) { // indices are used
            SourceVal = (double *)calloc(Vnum, sizeof(double));
        }
        
        for (int k = 0; k < Snum; k++) {
 
			if (mxGetM(prhs[1]) == 1) {
                int n = ((double *)mxGetData(prhs[1]))[k] - 1;
                for (i = 0; i < Vnum; i++) SourceVal[i] = mxGetInf();
                SourceVal[n] = 0.0;
            }
            else {
                SourceVal = (double *)mxGetData(prhs[1]) + Vnum*k;
            }
            
            double *DistanceData = (double *)mxGetData(plhs[0]) + Vnum*k;
			fmm->March(SourceVal, DistanceData);

        }
		if (mxGetM(prhs[1]) == 1) { free(SourceVal); }
    }
    else if (!strcmpi(command, "deinit")) {
	    FMM *fmm = (FMM*)((unsigned int *)mxGetPr(prhs[0]))[0];
		if (fmm) delete fmm;
    }
    else {
        mexPrintf("Invalid command");
    }
    
	return;
}
