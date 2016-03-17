/***************************************************************************/
/* MarchAlg.c  Fast Marching method for triangulations				   	   */
/***************************************************************************/

#include "y_stand.h" /* declaration file */
#include "mex.h"

static int *xStart,*yStart;							 /* Arrays that define the points that get zero intial 
														value, either all of them togheter or only one of
														them at a time.									*/
static int number_of_points;

static struct Vertex *V;							 /* The array of vertexes that is used in
													    the current fastmarch.							*/
static int Vnum;									 /* The number of vertexes in the array
													    of vertexes V.									*/


void insert(double v, int m);
void upheap(int k);
void remove_top(void);
int *BP;

/***************************************************************************/
/* SetSourceSet procedure, gets pointers to the arrays which are pointed   */
/* by xStart and yStart. npoints is the number of items in each array.	   */
/***************************************************************************/

void SetSourceSet(int xCoordinates[], int yCoordinates[], int points_no, struct Vertex VertexesArray[], int Vnum){
	static int i, ind;

	xStart = xCoordinates;
	yStart = yCoordinates;
	number_of_points = points_no;
	//matrix_allocate(sources_dis, Vnum, Vnum, double);


	for(i = 0; i < Vnum; i++){					/* Do some preperation work.			*/
		VertexesArray[i].source = FALSE;		/* Later the sources will be turned on. */
	}

											
	for (i = number_of_points - 1; i >=0; i--){ /* Turn on the the source indicator of the source variables.  */
		//ind = xStart[i]*YSIZE+yStart[i];
		ind = xStart[i];
		VertexesArray[ind].source = TRUE;
		//VertexesArray[ind].source_num = i;		/* Set the source num field.			*/
	}
	

}


/***************************************************************************/
/* IntiMarch procedure, Init phase of the procedure, the march is from the */
/* points specified by xplace and yplace and the index start_from that	   */
/* the algoritem gets. If the flag ALL_SOURCES_TOGETHER is defined the	   */
/* the procedure ignore the third argument and use all points in xplace	   */
/* and yplace as sources.												   */
/***************************************************************************/
static void
InitMarch(double srcval[], int Vnum){
	static int i;
	for(i = 0; i < Vnum; i++){

		if (srcval[i] >= MAX_DISTANCE) {
			V[i].U = MAX_DISTANCE;					/* Initialize the distance				*/
			BP[i] = Far;							/* Originally all are far.				*/
			V[i].source = FALSE;
			//V[i].source_num = 0;
		}
		else {
			V[i].U = srcval[i];
			insert(V[i].U, i);
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
static __forceinline boolean __fastcall
Quadratic(double qa, double qb, double qc, double *x1)
{

	static double d;

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
static __forceinline double __fastcall
update(struct Stencil *current_stencil, int k, int l){

	static double u1, u2;
	static double u;
	static double t;
	static double /*ctheta_mul_mb_div_ma_minus_one_mul_2,*/ ctheta_mul_ma_div_mb_minus_one_mul_2;
																			/* Both of them are between 0 and 1 becuase
																			   the triangles are not obuse.				*/
//	static double shortcut1;			/* Replaces 1 - a2 * (ctheta_mul_mb_div_ma_minus_one_mul_2 + 1) / b2 */
//	static double shortcut2;			/* Replaces Stheta * a2 */

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
TUpdateNeighbors(
	int    i,						/* The vertex to update				   */
	int	   becomes_alive			/* The vertex that now becomes alive.  */
){
	double          u, t3, u0;
	int             n = BP[i];
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
			if(BP[k] == Alive || BP[l] == Alive){		/* Do update only if the two other vertexes
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
		insert(u, i);	/* also changes BP to Trail	   */
		V[i].U = u;
	}
	else {				/* change position in the heap */
		if (a[n].u > u) {
			a[n].u = V[i].U = u;
			upheap(n);
		}
	}
}
/***************************************************************************/
/* Update Triangle Neighbors procedure					   */
/***************************************************************************/
static void __forceinline
__fastcall CallTUpdateNeighbors(
	int             i,	/* vertex number */
	int start_from		/* The source from which the current march began.	*/
){
	static int local_vn_count;			/* The number of neighbours yet to cover. */
	static int *VN_ind;					/* Pointer to the array of neighbours */
										/* Call TUpdateNeighbors for every neighbour vertex of
										   the given vertex. */

//	mexPrintf("%d  (%d)\n", i,start_from);

	for (VN_ind = V[i].VN, local_vn_count = 0; local_vn_count < V[i].vn; VN_ind++, local_vn_count++)
										/* If the neighbour is not alive					*/
		if(BP[*VN_ind] is_not Alive /*&& V[*VN_ind].source == FALSE*/){
			TUpdateNeighbors(*VN_ind, i);
		}

}

/***************************************************************************/
/* MarchAlg procedure, for fast marching method, extentions of Sethian     */
/* This procedure works with triangulations that contains only non-obtuse  */
/* triangles															   */
/***************************************************************************/
void __fastcall MarchAlg(
	struct Vertex  *temp_V,	/* vertices coordinates 1..Vnum */
	int             temp_Vnum,
	double *srcval){
	V = temp_V;						/* Save the array of vertexes in the V array.		*/
	Vnum = temp_Vnum;				/* Save the number of vertexes in the array V.		*/

	InitMarch(srcval, Vnum);

									/* CAllTUpdateNeighbours for every vertex when it becomes
									   the smallest vertex in the heap. Exit the loop when
									   there are no more sources which are not alive */
	


	while (N != 0){
		CallTUpdateNeighbors(a[1].v, 0);
		remove_top();
	}

#ifndef NO_MANIFOLD
	{
		int             k;
		for (k = temp_Vnum - 1; k >= 0; temp_V++, k--) /* Do for every vertex */
			if (temp_V->U > MAX_DISTANCE / 2)
				temp_V->U = 0;
	}
#endif

}

