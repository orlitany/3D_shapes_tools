/***************************************************************************/
/* my_heap.c : Fast marching method - The Eikonal triangular case 	       */
/* This program was written by Ronny Kimmel    							   */
/***************************************************************************/

#include "y_stand.h"    /* declaration file */

/***************************************************************************/
/*  Heap := Parent of j is  at j div 2, children of j are at 2j and 2j+1   */
/* Min heap: Parent is always less than its children                       */
/* the root is at a[1].u, a[0].u is an anchor value                        */
/***************************************************************************/


/***************************************************************************/
/* initialize the heap to be an empty heap with an anchor value only.	   */
/***************************************************************************/
void Heap::init_heap(int Vnum){
	/* The heap size is at most the number of vertices Vnum, allocate the
	   heap array and the back pointer array.							   */
    a = (struct heap_element *) malloc(sizeof(struct heap_element) *(Vnum+1));
    BP = (int *) malloc(sizeof(int) *Vnum);
	
	a[0].u = -MAX_DISTANCE;			/* anchor the root of the heap		   */
	N=0;							/* Originaly the heap is empty.		   */
}

/***************************************************************************/
/* free the memory used by the heap.									   */
/***************************************************************************/
void Heap::free_heap(void){
	free(a);
	free(BP);
}

void Heap::reset(void) {
	a[0].u = -MAX_DISTANCE;			/* anchor the root of the heap		   */
	N=0;							/* Originaly the heap is empty.		   */
}


/***************************************************************************/
/* upheap - gets the number of a vertex and moves it up the heap structure */
/* till it is bigger than its father.									   */
/***************************************************************************/
void
Heap::upheap(int k)
{
	double          u = a[k].u;	/* current value */
	int             v = a[k].v;

	while (a[k >> 1].u > u) {	/* While the place was not found move the vertex up */
		a[k] = a[k >> 1];
		BP[a[k].v] = k;	        /* back pointer */
		k >>= 1;
	}
	a[k].u = u;
	a[k].v = v;
	BP[a[k].v] = k;	/* back pointer */
}
/***************************************************************************/
/* insert an elemet to the heap structure								   */
/* v is the number of the vertex to add in the V array and u is its U value*/
/***************************************************************************/
void
Heap::insert(double u, int v)
{
	N++;
	a[N].u = u;
	a[N].v = v;
	BP[v] = N;		/* back pointer */
	upheap(N);
}

/***************************************************************************/
/* downheap - get the number of a vertex in the heap and moves it down the */
/* heap structure till its two sons are bigger than him.				   */
/***************************************************************************/
void
Heap::downheap(int k)
	/* k index to the element to be moved down the
						   heap */
{
	double          u = a[k].u;	/* current value */
	int             v = a[k].v;
	int             j = k;				/* Moves down the heap, till a place for vertex v is found. */

	while ((j <<= 1) <= N) {			/* While have not reached the button and cont is still true */
		if (j < N and a[j].u > a[j + 1].u) /* If the brother of j is smaller, move j to the brother.*/
			j++;
		if (u <= a[j].u)				/* If vertex k is smaller than vertex j end search.			*/
			break;
		else {
			a[k] = a[j];				/* Move vertex k one plece ower.							*/
			BP[a[k].v] = k;				/* update back pointer */
			k = j;		
		}
	}
	a[k].u = u;							/* We have found a place for v, place the vertex here.		*/
	a[k].v = v;
	BP[a[k].v] = k;	/* back pointer */
}
/***************************************************************************/
/* remove_top - removes the vertex with the minimum U value from the heap  */
/* structure															   */
/***************************************************************************/
void
Heap::remove_top()
{
	BP[a[1].v] = Alive;	/* set back pointer to Alive (not in list)  */
	a[1] = a[N];		/* Move the last vertex to be the first     */
	BP[a[1].v] = 1;		/* back pointer								*/
	N--;
	downheap(1);		/* Update the heap structure.				*/
}
