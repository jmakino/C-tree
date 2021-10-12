// 
// BHtree.C
// Version 1999/1/1 --- cleaned up and some comments  added.
//
// The tree construction/handling package designed for TREE
// implementation of SPH/NBODY programs
//
//
// This program uses a new tree construction method, based on
// Morton ordering and top-down tree construction.
//
// non-local functions (not complete yet...)
//
// setup_tree()
//    allocate the memory for tree if necessary
//    create the tree structure
// set_cm_quantities_for_default_tree()
//    calculate mass and position for tree nodes
//
// set_hmax_for_sph() (SPH use only)
// check_and_set_nbl() (SPH use only)
//
// calculate_gravity_using_tree() (NO GRAPE)
//
// evaluate_gravity_using_tree_and_list()(GRAPE)

extern "C" double cpusec();

#ifndef PRC
#define PR(x)  cerr << #x << " = " << x << " "
#define PRC(x) cerr << #x << " = " << x << ",  "
#define PRL(x) cerr << #x << " = " << x << "\n"
#endif

#include  <stdlib.h>
#include  <math.h>
#include  <iostream>
#include "utils.hpp"

const int LISTMAX = 40000;
using namespace std;


#define real double
#include "BHtree.h"
#include "nbody_system.h"

int bhnode::nplimit = 1;

//void  nbody_particle::correct_phi_self_gravity(real epsinv)
//{phi_gravity += mass*epsinv;}
//void  nbody_particle::clear_acc_phi_gravity()
//{acc_gravity = 0.0;phi_gravity = 0.0;}

void dump_octal(BHlong x)
{
    char st[256];
    sprintf(st," %llo ",x);
    cerr <<  st ;
}

    

BHlong conv_to_morton_old(int ix, int keybits)
{
    BHlong dum = 0;
    //    char st[256];
    //    cerr << "conv_to_morton "; PR(ix); PRL(keybits);
    //    sprintf(st,"%lo",ix);
    //    cerr << "ix = " << st << endl;
    int i, j;
    for(i = j= 0; i<keybits; i++,j+=3){
	if (ix & (1<<i)){
	    dum |= ((BHlong) 1)<<j;
	}
    }
    //sprintf(st,"%lo",dum);
    //    cerr << "dum = " << st << endl;
    return dum;

}
BHlong conv_to_morton(int ix, int keybits)
{
    BHlong dum = 0;
    //    char st[256];
    //    cerr << "conv_to_morton "; PR(ix); PRL(keybits);
    //    sprintf(st,"%lo",ix);
    //    cerr << "ix = " << st << endl;
    int i, j;
    for(i = j= 0; i<keybits; i++,j+=3){
	if ((ix>>i) & 1){
	    dum |= ((BHlong) 1)<<j;
	}
    }
    //sprintf(st,"%lo",dum);
    //    cerr << "dum = " << st << endl;
    return dum;

}
    
inline BHlong  construct_key(const vector & pos, real rscale, int ioffset, int keybits)
{
    BHlong ix[3];
    vector p = pos;
    for(int i = 0; i<3; i++){
	ix[i] = (BHlong) (p[i]*rscale+ioffset);
    }
    return (conv_to_morton(ix[0],keybits)<<2 )
	|(conv_to_morton(ix[1],keybits)<<1 )
	|(conv_to_morton(ix[2],keybits));
}
    

void bhparticle::set_key(real rscale, int ioffset, int keybits)
{
    key =  construct_key(rp->get_pos(), rscale, ioffset, keybits);
    //    PRL(key);
    
}

int compare_key(const void * p1, const void * p2)
{
    bhparticle * pp= (bhparticle*)p1;
    bhparticle * qq=(bhparticle*)p2;
    BHlong comp = ((BHlong) pp->get_key()) - ((BHlong) qq->get_key());
    if (comp > 0L){
	return 1;
    }else if (comp == 0L){
	return 0;
    }else{
	return -1;
    }
}

void sort_bh_array( bhparticle * r, int lo, int up )
{
    int i, j;
    bhparticle tempr;
    while ( up>lo ) {
	i = lo;
	j = up;
	tempr = r[lo];
	/*** Split file in two ***/
	while ( i<j ) {
	    for ( ; r[j].key > tempr.key; j-- );
	    for ( r[i]=r[j]; i<j && r[i].key<=tempr.key; i++ );
	    r[j] = r[i];
	}
	r[i] = tempr;
	/*** Sort recursively, the smallest first ***/
	if ( i-lo < up-i ) { sort_bh_array(r,lo,i-1);  lo = i+1; }
	else    { sort_bh_array(r,i+1,up);  up = i-1; }
    }
}
void check_bh_array( bhparticle * r, int size )
{
    for(int i = 0; i<size-1;i++){
	if(r[i].get_key() > r[i+1].get_key()){
	    PR(i); PR(r[i].get_key()); PRL(r[i+1].get_key());
	    cerr << "Sort failed ... \n";
	    exit (1);
	}
    }
}

real initialize_key(int nbody,
		    real_particle * rp,
		    int & nkeysize,
		    bhparticle * &bhp)
{
    if (nbody > nkeysize || bhp == NULL){
	if (bhp != NULL){
	    delete [] bhp;
	}
	nkeysize = nbody+100;
	bhp = new bhparticle[nkeysize];
#ifdef REUSE_PREVIOS_DATA
	// With present quicksort routine, the pre-sorted data
	// actually DEGRADE its performance. So DO NOT ACTIVATE
	// THIS PART --- JM 1998/12/22
	for(int i = 0; i<nbody; i++){
	    bhparticle * p = bhp + i;
	    p->set_rp(rp+i);
	}
#endif	
    }
    real rmax = 1;
    for(int i = 0; i<nbody; i++){
	vector p = (rp+i)->get_pos();
	for (int k = 0; k<3; k++){
	    if (fabs(p[k])>=rmax) rmax *= 2;
	}
    }
    real rscale = 1.0/rmax*default_ix_offset;
	
    for(int i = 0; i<nbody; i++){
	bhparticle * p = bhp + i;
#ifndef REUSE_PREVIOS_DATA	
	p->set_rp(rp+i);
#endif	
	p->set_key(rscale, default_ix_offset, default_key_length);
	//	PR(i); PRL(p->get_key());
    }
    cerr << "Call quicksort, cpu = " <<NbodyLib::GetWtime() << endl;
    //    sort_bh_array(bhp,0,nbody-1);
    qsort(bhp, nbody, sizeof(bhparticle), compare_key);
    // The private sort routine is for some unknow reason
    // much faster than qsort of the system for large N
    cerr << "Exit quicksort, cpu = " <<NbodyLib::GetWtime() << endl;
    for(int i = 0; i<nbody; i++){
	bhparticle * p = bhp + i;
	//	PR(i); PR(p->get_key()); PRL(p->get_rp()->get_index());
    }
    return rmax;
}

void bhnode::assign_root(vector root_pos, real length, bhparticle * bp, int np)
{
    pos = root_pos;
    l = length;
    bpfirst = bp;
    nparticle = np;
}


    


void bhnode::create_tree_recursive(bhnode * & heap_top, int & heap_remainder,
				   BHlong current_key,
				   int current_level,
				   int n_critical)
{
//    cerr << "create tree called "; PRC(nparticle); PRL(n_critical);
//    PRL(heap_remainder);
    if (heap_remainder <= 0){
	cerr << "create_tree: no more free node... exit\n";
	exit(1);
    }
    if (nparticle <= n_critical) return;
    if (current_level == 0) return;
    //    cerr << "Enter recursion\n";
    //    dump();
    BHlong keyscale = ((BHlong) 1)<<((current_level-1)*3);
    bhparticle * bptmp = bpfirst;
    int npremain = nparticle;
    for(int i=0; i<8;i++)child[i] = NULL;
    isleaf = 1;
    for(int i=0; i<8;i++){
	BHlong new_key = current_key + keyscale * i;
	vector new_pos = pos + vector( ((i&4)*0.5-1)*l/4,
				       ((i&2)    -1)*l/4,
				       ((i&1)*2  -1)*l/4);
	
	if(bptmp->get_key() - new_key <keyscale){
	    // current bptmp is actually in the current subnode
	    // search for the end location
	    int p0 = 0;
	    int p1 = npremain-1;
	    if ((bptmp+p1)->get_key() - new_key >=keyscale){
		while (p1 - p0 > 1){
		    int pnew = (p0+p1)/2;
		    if ((bptmp+pnew)->get_key() - new_key <keyscale){
			p0 = pnew;
		    }else{
			p1 = pnew;
		    }
		}
		p1 = p0;
	    }
	    p1 ++;
	    isleaf = 0;
	    child[i] = heap_top;
	    heap_top ++;
	    heap_remainder -- ;
	    child[i]->bpfirst = bptmp;
	    child[i]->pos = new_pos;
	    child[i]->l = l*0.5;
	    child[i]->nparticle = p1;
	    child[i]->isleaf = 1;
	    child[i]->create_tree_recursive(heap_top, heap_remainder,
					    new_key, current_level-1, n_critical);
	    bptmp += p1;
	    npremain -= p1;
	    if (npremain <= 0) return;
				  
	}
    }

    //dump();
}


void spc(int indent)
{
    for(int i=0;i<indent;i++)printf(" ");
}

void bhnode::dump(int indent)
{
    int i;
    spc(indent); printf("node pos %g %g %g", pos[0], pos[1], pos[2]);
    printf(" cmpos %g %g %g", cmpos[0], cmpos[1], cmpos[2]);
    printf(" mass %g", cmmass);

    if (isleaf){
	spc(indent); printf(" LEAF np=%d\n", nparticle);
	bhparticle * bp = bpfirst;
	for(i = 0; i < nparticle; i++){
	    spc(indent+2);
	    real_particle * p = (bp+i)->get_rp();
	    vector x=p->get_pos();
	    printf("index=%d, pos=%g %g %g\n", p->get_index(),
		  x[0], x[1], x[2]);
	}
    }else{
	spc(indent); printf(" NON LEAF\n");
	for(i=0;i<8;i++){
	    if (child[i] != NULL){
		child[i]->dump(indent + 2);
	    }
	}
    }
}

// inbox: returns 0 if the particle is in the box;

int inbox(vector  & cpos, // center of the box
	  vector  & pos,  // position of the particle
	  real l)         // length of one side of the box
    
{
    for(int  i = 0; i< ndim; i++){
	if (fabs(pos[i]-cpos[i]) > l*0.5) return 1;
    }
    return 0;
}
	
	
int bhnode::sanity_check()
{
    int i;
    int iret = 0;
    if (isleaf){
	// this is the lowest level node. Things to check:
	// all particles are in the cell
	bhparticle * bp = bpfirst;
	cout << "Leaf np="<< nparticle <<endl;
	for(i = 0; i < nparticle; i++){
	    real_particle * p = (bp+i)->get_rp();
	    vector ppos = p->get_pos();
	    if(inbox(pos,ppos,l)){
		cerr << "Error, particle out of box ... \n";
		dump();
		return 1;
	    }
	}
    }else{

	// This is the non-leaf node. Check the position and side
	// length of the child cells and then check recursively..
	cout << "Non Leaf " << pos  <<endl;
	for(i=0;i<8;i++){
	    if (child[i] != NULL){
		int err = 0;
	        err = child[i]->sanity_check();
		if (l*0.5 != child[i]->get_length()) err += 2;
		vector relpos = pos-child[i]->get_pos();
		for (int k = 0 ; k<ndim;k++){
		    if (fabs(relpos[k]) !=l*0.25)err += 4;
		}
		if (err){
		    cerr << "Child " << i << " Error type = " << err << endl;
		    dump();
		}
		iret += err;
	    }
	}
    }
    return iret;
}

void bhnode::setup_child_nodes(child_nodes* &c,
			       child_nodes* cbase,
			       bhnode* nfirst,
			       bhparticle* pfirst)
{
    int i;
    //    printf("setup at address %d\n", (int)(c-cbase));
    auto myc = c;
    if (isleaf){
	// this is the lowest level node. Infomation is already set at upper
	// level
	return;
    }else{
	// This is the non-leaf node. Check the position and side
	// length of the child cells and then check recursively..
	//	cout << "Non Leaf " << pos  <<endl;
	for(i=0;i<NCHILDLEN;i++){
	    c->isleaf[i] = 1;
	    c->nparticle[i] = 0;
	}
	    
	for(i=0;i<NCHILDLEN;i++){
	    if (child[i] != NULL){
		c->posx[i] = child[i]->pos[0];
		c->posy[i] = child[i]->pos[1];
		c->posz[i] = child[i]->pos[2];
		c->l[i] = child[i]->l;
		c->cmposx[i] = child[i]->cmpos[0];
		c->cmposy[i] = child[i]->cmpos[1];
		c->cmposz[i] = child[i]->cmpos[2];
		c->cmmass[i] = child[i]->cmmass;
		c->cindex[i] = -1;
		c->pindex[i] = child[i]->bpfirst -pfirst;
		c->nparticle[i] = child[i]->nparticle;
		c->isleaf[i] = child[i]->isleaf;
		// if (c->isleaf[i] ){
		//     printf("Leaf %d %g %g %g %g %d\n",
		// 	   i, c->posx[i], c->posy[i], c->posx[i],c->cmmass[i],
		// 	   c->nparticle[i]);
		// }
			   
	    }else{
		c->isleaf[i] = 1;
		c->nparticle[i] = 0;
	    }
	}
	for(i=0;i<NCHILDLEN;i++){
	    if (child[i] != NULL && !child[i]->isleaf){
		c+=1;
		myc->cindex[i] = c-cbase;
		child[i]->setup_child_nodes(c,cbase, nfirst,pfirst);
	    }
	}
	// auto myid = myc - cbase;
	// for(i=0;i<NCHILDLEN;i++){
	//     if (myc->cindex[i]>0){
	// 	printf("my=%d cindex[%d]=%d\n", (int)myid,i, (int)(myc->cindex[i]));
	//     }
	//     if (myc->nparticle[i]){
	// 	printf("my=%d nparticle[%d]=%d\n", (int)myid,i, (int)(myc->nparticle[i]));
	//     }
	// }
    }
    c++;
}
void nbody_system::dump_child_nodes(child_nodes * cbase,
				    int ic,
				    bhnode* nfirst,
				    bhparticle* pfirst,
				    int level)
{
    int i;
    child_nodes * c = cbase + ic;
    for(i=0;i<NCHILDLEN;i++){
	if (!c->isleaf[i] || c->nparticle[i]){
	    for(auto is=0;is<level; is++)printf(" ");
	    printf("ic:%d i:%d cpos: %g %g %g cm: %g %g %g m: %g l:%g id:%ld np:%ld\n",
		   ic, i, c->posx[i], c->posy[i], c->posz[i], 
		   c->cmposx[i], c->cmposy[i], c->cmposz[i],
		   c->cmmass[i], c->l[i], c->cindex[i], c->nparticle[i]);
	    if (!c->isleaf[i]){
		dump_child_nodes(cbase, c->cindex[i], nfirst, pfirst, level+2);
	    }
	}
	if (c->isleaf[i] && c->nparticle[i]){
	    for (auto ip=0;ip<c->nparticle[i];ip++){
		for(auto is=0;is<level+2; is++)printf(" ");
		real_particle * p = (pfirst + c->pindex[i]+ip)->rp;
		printf("id:%d ppos: %g %g %g  m: %g\n",
		       p->get_index(), p->pos[0],  p->pos[1],   p->pos[2], p->mass);
	    }
	}
    }
}

void nbody_system::printchild_nodes(child_nodes * cbase,
		      int nc)
{
    int i;
    printf("Print child called for %d npdes\n", nc);
    for(i=0;i<nc;i++){
	int sum=0;
	auto c = cbase+i;
	for(auto k=0;k<NCHILDLEN; k++)sum+= c->nparticle[k];
	for(auto k=0;k<NCHILDLEN; k++){
	    if (c->isleaf[k]==0 || c->nparticle[k]){
		
		printf("i:%d cpos: %g %g %g cm: %g %g %g m: %g l: %g %ld %ld %ld\n",
		       i,c->posx[k], c->posy[k], c->posz[k], 
		       c->cmposx[k], c->cmposy[k], c->cmposz[k],
		       c->cmmass[k], c->l[k], c->cindex[k],
		       c->pindex[k], c->nparticle[k]);
	    }
	}
    }
}

#ifdef SPH	
void  bhnode::set_hmax_for_sph()
{
    int i;
    hmax_for_sph = 0;
    if (isleaf){
	bhparticle * bp = bpfirst;
	for(i = 0; i < nparticle; i++){
	    real hp = (bp+i)->get_rp()->get_h();
	    if(hmax_for_sph < hp) hmax_for_sph = hp;
	}
    }else{
	for(i=0;i<8;i++){
	    if (child[i] != NULL){
		child[i]->set_hmax_for_sph();
		if (hmax_for_sph < child[i]->hmax_for_sph)
		    hmax_for_sph = child[i]->hmax_for_sph;
	    }
	}
    }
}
#endif

void  bhnode::set_cm_quantities()
{
    int i;
    cmpos = 0.0;
    cmmass = 0.0;
    if (isleaf){
	bhparticle * bp = bpfirst;
	for(i = 0; i < nparticle; i++){
	    real mchild = (bp+i)->get_rp()->get_mass();
	    cmpos += mchild*(bp+i)->get_rp()->get_pos();
	    //	    vector  p = (bp+i)->get_rp()->get_pos();
	    //	    cmpos[0] += mchild*p[0];
	    //	    cmpos[1] += mchild*p[1];
	    //	    cmpos[2] += mchild*p[2];
	    cmmass += mchild;
	}
    }else{
	for(i=0;i<8;i++){
	    if (child[i] != NULL){
		child[i]->set_cm_quantities();
		real mchild = child[i]->cmmass;
		cmpos += mchild*child[i]->cmpos;
		cmmass += mchild;
	    }
	}
    }
    cmpos /= cmmass;
}

inline real separation_squared(bhnode * p1,bhnode * p2)
{
    real r2 = 0;
    real xmin = (p1->get_length()+ p2->get_length())*0.5;
    vector dx = p1->get_pos() - p2->get_pos();
    for (int k = 0; k<ndim; k++){
	real adx = fabs(dx[k]);
	if (adx > xmin){
	    adx -= xmin;
	}else{
	    adx = 0;
	}
	r2 += adx*adx;
    }
    return r2;
}

inline real separation_squared0(bhnode * p1, vector & pos2)
{
    real r2 = 0;
    real xmin = p1->get_length()*0.5;
#if 0
    vector dx = p1->get_pos() - pos2;
    for (int k = 0; k<ndim; k++){
	real adx = fabs(dx[k]);
	if (adx > xmin){
	    adx -= xmin;
	    r2 += adx*adx;
	}
    }
#else
    real adx, ady, adz;
    vector * p = p1->get_posp();
    adx = fabs((*p)[0]-pos2[0]);
    ady = fabs((*p)[1]-pos2[1]);
    adz = fabs((*p)[2]-pos2[2]);
    if (adx>xmin){
	adx -= xmin;
	r2 += adx*adx;
    }
    if (ady>xmin){
	ady -= xmin;
	r2 += ady*ady;
    }
    if (adz>xmin){
	adz -= xmin;
	r2 += adz*adz;
    }
#endif	       
    
    return r2;
}

inline real separation_squared_old(bhnode * p1, vector & pos2)
{
    real r2 = 0;
    real xmin = p1->get_length()*0.5;
    real adx, ady, adz;
    vector * p = p1->get_posp();
    adx = fabs((*p)[0]-pos2[0])-xmin;
    ady = fabs((*p)[1]-pos2[1])-xmin;
    adz = fabs((*p)[2]-pos2[2])-xmin;
    if (adx>0.0){
	r2 += adx*adx;
    }
    if (ady>0.0){
	r2 += ady*ady;
    }
    if (adz>0.0){
	r2 += adz*adz;
    }
    
    return r2;
}
inline real separation_squared(bhnode * p1, vector & pos2)
{
    real r2 = 0;
    auto dx = p1->get_pos()-pos2;
    return dx*dx;
}

inline int  are_overlapped(bhnode * p1,bhnode * p2)
{
    real xmin = (p1->get_length()+ p2->get_length())*0.4999999999999999;
    vector dx = p1->get_pos() - p2->get_pos();
#if 0
    for (int k = 0; k<ndim; k++){
	if(fabs(dx[k]) > xmin) return 0;
    }
#else
    if((fabs(dx[0]) > xmin) ||(fabs(dx[1]) > xmin) ||(fabs(dx[2]) > xmin)) {
	return 0;
    }
	
#endif
    return 1;
}


inline int  are_overlapped_with_cutoff(bhnode * p1,bhnode * p2,
				       real cutoff)
{
    real xmin = (p1->get_length()+ p2->get_length())*0.4999999999999999
	+cutoff;
    vector dx = p1->get_pos() - p2->get_pos();
    if((fabs(dx[0]) > xmin) ||(fabs(dx[1]) > xmin) ||(fabs(dx[2]) > xmin)) {
	return 0;
    }
    real xmin2 = (p1->get_length()+ p2->get_length())*0.4999999999999999;
    real r2sum = 0;
    for(int k=0;k<3;k++){
	real d = fabs(dx[k]) - xmin2;
        if (d > 0) r2sum += d*d;
    }
    if (r2sum > cutoff*cutoff) {
	return 0;
    }
    return 1;
}
inline int  are_overlapped_with_cutoff(bhnode * p1, bhparticle * p2,
				       real cutoff)
{
    real xmin = p1->get_length()*0.4999999999999999+cutoff;
    vector dx = p1->get_pos() - p2->get_rp()->get_pos();
    if((fabs(dx[0]) > xmin) ||(fabs(dx[1]) > xmin) ||(fabs(dx[2]) > xmin)) {
	//	cerr << "returns 0\n";
	return 0;
    }
    //    PRC(p1->get_pos());PRC(p1->get_length());
    //    PRC(p2->get_rp()->get_pos());
    real xmin2 = p1->get_length()*0.4999999999999999;
    real r2sum = 0;
    for(int k=0;k<3;k++){
	real d = fabs(dx[k]) - xmin2;
        if (d > 0) r2sum += d*d;
    }
    if (r2sum > cutoff*cutoff) {
	//	cerr << "returns 0\n";
	return 0;
    }
    //    cerr << "returns 1\n";
    return 1;
}

#ifdef SPH
int check_and_set_nbl(bhnode * p1,bhnode * p2);
int check_and_set_nbl(bhnode * p1,bhnode * p2)
{
    int iret = 0;
    if(p1 == NULL) return iret;
    if(p2 == NULL) return iret;
    real rcrit = p1->hmax_for_sph;
    if (rcrit <  p2->hmax_for_sph)rcrit = p2->hmax_for_sph;
    rcrit *= 2;
    real rcrit2 = rcrit*rcrit;
    if(separation_squared(p1,p2) > rcrit2) return iret;
    if (p1->isleaf == 0 || p2->isleaf == 0){
	//
	// either node is not leaf. Go down
	//
	if (p1->isleaf || (p2->isleaf == 0) && p2->l > p1->l){
	    //
	    // p1 is already leaf, or both are node but p2 is larger.
	    // go down p2 by 
	    for(int i = 0; i<8;i++)
		iret |= check_and_set_nbl(p1, p2->child[i]);
	}else{
	    //
	    // now, we can go down p1 ...
	    for(int i = 0; i<8;i++)
		iret |= check_and_set_nbl(p1->child[i], p2);
	}
    }else{
	//
	// both are leaves; can directly evaluate particles
	bhparticle * bpi = p1->bpfirst;
	bhparticle * bpj = p2->bpfirst;
	for(int i = 0; i < p1->nparticle; i++)
	    for(int j = 0; j < p2->nparticle; j++){
		if((bpi+i)->get_key() <(bpj+j)->get_key() ){
		    //
		    // here, comparison of key guarantee that a pair
		    // is evaluated only once
	  
		    iret |=check_and_set_nbl(*((bpi+i)->get_rp()),
					     *((bpj+j)->get_rp()));
		}
	    }
    }
    return iret;
}

#endif
#if 0
static bhparticle * bp = NULL;
static int bhpsize = 0;
static int bnsize = 0;
static bhnode * bn;
#endif
void real_system::set_cm_quantities_for_default_tree()
{
    bn->set_cm_quantities();
}


void real_system::setup_tree()
{
    real rsize = initialize_key(n,get_particle_pointer(),bhpsize,bp);
    cout << "Setup tree: called\n";
    int expected_bnsize =  (int)(bhpsize*0.6+100);
    if (bnsize < expected_bnsize){
	if (bnsize != 0){
	    delete [] bn;
	    delete [] cn;
	}
	bnsize = expected_bnsize;
	bn = new bhnode[bnsize];
	//	cn = new child_nodes[bnsize/4+100];
	cn = new child_nodes[bnsize];
    }
    for(int j = 0; j<bnsize;j++) (bn+j)->clear();
    bn->assign_root(vector(0.0), rsize*2, bp, n);
    bhnode * btmp = bn+1;
    int heap_remainder = bnsize-1;
    BHlong key = 0;
    bn->create_tree_recursive(btmp,heap_remainder,key,
			      default_key_length, 8 );
    PR(bnsize);    PRL(heap_remainder);
    //    PRL(bn->sanity_check());
}

	
#ifdef SPH	
int sph_system::set_nnb_using_tree()
{
    setup_tree();
    real_particle * psph = get_particle_pointer();
    apply_vf(real_particle::clear_nnb);
    bn->set_hmax_for_sph();
    int iret = check_and_set_nbl(bn, bn);
    apply_vf(real_particle::sort_nblist);
    return iret;
}
#endif

inline void accumulate_force_from_point(vector dx, real r2, real eps2, 
				 vector & acc,
				 real & phi,
				 real jmass)
{
    double r2inv = 1/(r2+eps2);
    double rinv  = sqrt(r2inv);
    double r3inv = r2inv*rinv;
    phi -= jmass*rinv;
    acc += jmass*r3inv*dx;
}

static real total_interactions;
static int tree_walks;
static int nisum;
static double t_walk;
static double t_calc;
void clear_tree_counters()
{
    total_interactions = 0;
    tree_walks = 0;
    nisum = 0;
    t_walk = t_calc = 0;
}
void print_tree_counters()
{
    real avg = total_interactions/nisum;
    PRC(nisum); PRC(tree_walks); PRC(total_interactions); PRL(avg);
    cout <<"tree_walks = " <<tree_walks << " ntaverage = " << avg << endl;
    cout <<"walk time = " <<t_walk << " calc time = " << t_calc << endl;
}

void calculate_force_from_interaction_list(const vector & pos,
					   real eps2, 
					    vector & acc,
					    real & phi,
					    vector * poslist,
					    real * masslist,
					    int list_length)
{
    acc = 0.0;
    phi = 0.0;
    for(int i = 0; i<list_length; i++){
	vector dx = *(poslist+i)-pos;
	real r2 = dx*dx;
	accumulate_force_from_point(dx, r2, eps2, acc, phi,*(masslist+i));
    }
}



void bhnode::accumulate_force_from_tree(vector & ipos, real eps2, real theta2,
					vector & acc,
					real & phi)
{
    vector dx = cmpos - ipos;
    real r2 = dx*dx;
    if (r2*theta2 > l*l){
	// node and position is well separated;
	accumulate_force_from_point(dx, r2, eps2, acc, phi, cmmass);
	total_interactions += 1;
    }else{
	int i;
	if (isleaf){
	    bhparticle * bp = bpfirst;
	    for(i = 0; i < nparticle; i++){
		vector dx = (bp+i)->get_rp()->get_pos()-ipos;
		real r2 = dx*dx;
		accumulate_force_from_point(dx, r2, eps2, acc, phi,
					    (bp+i)->get_rp()->get_mass());
		total_interactions += 1;
	    }
	}else{
	    for(i=0;i<8;i++){
		if (child[i] != NULL){
		child[i]->accumulate_force_from_tree(ipos,eps2,theta2, acc, phi);
		}
	    }
	}
    }
}

void bhnode::add_to_interaction_list(bhnode & dest_node, real theta2,
				     vector * pos_list,
				     real * mass_list,
				     int & nlist,
				     int list_max,
				     int & first_leaf)
{
#if 0
    printf("add: nlist=%d pos, l= %g %g %g %g\n",
	   nlist,  pos[0],pos[1],pos[2], l );
    printf("r2*theta2= %g l*l=%g\n",
	   separation_squared(&dest_node,cmpos),
	   l*l);
#endif    
    if((separation_squared(&dest_node,cmpos)*theta2 > l*l)
       && (!are_overlapped(this,&dest_node) ) ){
	// node and position is well separated;
	//	printf("opened\n");
	*(pos_list+nlist) = cmpos;
	*(mass_list+nlist) = cmmass;
	nlist ++;
	if (nlist > list_max){
	    cerr << "List length exceeded\n";
	    exit(1);
	}
	
    }else{
	int i;
	if (isleaf || (this == (&dest_node))){
	    if (this == (&dest_node)){
		// adding the particles in the node itself
		first_leaf = nlist;
	    }
	    bhparticle * bp = bpfirst;
	    for(i = 0; i < nparticle; i++){
		*(pos_list+nlist) = (bp+i)->get_rp()->get_pos();
		
		*(mass_list+nlist) =(bp+i)->get_rp()->get_mass();
		nlist ++;
		if (nlist > list_max){
		    cerr << "List length exceeded\n";
		    exit(1);
		}
	    }
	}else{
	    for(i=0;i<8;i++){
		if (child[i] != NULL){
		    child[i]->add_to_interaction_list(dest_node, theta2,
						      pos_list, mass_list,
						      nlist, list_max,
						      first_leaf);
		}
	    }
	}
    }
}
template<class T>
inline auto absfloor(T x, T f)
{
    if (x < 0) x= -x;
    auto x1 = x - f;
    if (x1 > 0){
	return x1;
    }else{
	return (T) 0;
    }
}

template<class T>
inline int64_t absge_int64(T x, T f)
{
    if (x < 0) x= -x;
    if (x >= f){
	return (int64_t)1;
    }else{
	return (int64_t)0;
    }
}
template<class T>
inline bool absge(T x, T f)
{
    if (x < 0) x= -x;
    return (x >= f);
}

inline void cn_opening_criterion(child_nodes *c,
				 vector & posi,
				 real li,
				 real theta2,
				 int64_t well_separated[NCHILDLEN])
{
    //#pragma clang loop vectorize_width(2,scalable)
    //#pragma clang loop vectorize(enable)
#pragma omp simd
    for (auto j=0; j< NCHILDLEN; j++){
	auto  xmin = (c->l[j] + li)*0.499999999999999;
	auto xout = absge(c->posx[j] - posi[0], xmin);
	auto yout = absge(c->posy[j] - posi[1], xmin);
	auto zout = absge(c->posz[j] - posi[2], xmin);
	auto dx = c->cmposx[j] -posi[0];
	auto dy = c->cmposy[j] -posi[1];
	auto dz = c->cmposz[j] -posi[2];
	auto r2 = dx*dx + dy*dy+ dz*dz;
	if (r2* theta2 > c->l[j]*c->l[j] && (xout ||yout||zout)){
	    well_separated[j]=1;
	} else{
	    well_separated[j]=0;
	}
#if 0	
	printf("dx = %g %g %g %g r2*theta2= %g l*l=%g\n",
	       dx, dy, dz, r2, r2*theta2, c->l[j]*c->l[j] );
	printf("node %d-%d pos: %g %g %g m:%g l:%g wellsep=%d\n",
	       cnindex, j, c->posx[j], c->posy[j], c->posz[j],
	       c->cmmass[j], c->l[j], well_separated[j]);
#endif	
    }
}
void cn_add_to_interaction_list(child_nodes * cn,
				int cnindex,
				bhparticle * bn,
				bhnode & dest_node,

				real theta2,
				vector * pos_list,
				real * mass_list,
				int & nlist,
				int list_max,
				int & first_leaf)
{
    bool well_separated[NCHILDLEN];
    vector posi = dest_node.get_pos();
    real   li   = dest_node.get_length();
    bhparticle* b = dest_node.bpfirst;
    auto c = cn+cnindex;
#if 0    
    printf("cn_add: cnindex=%d nlist=%d pos, l= %g %g %g %g\n",
	   cnindex, nlist, posi[0],posi[1],posi[2], li );
#endif

    //    cn_opening_criterion(c,posi,li,theta2,well_separated);
#pragma clang loop vectorize(enable)
    for (auto j=0; j< NCHILDLEN; j++){
	auto  xmin = (c->l[j] + li)*0.499999999999999;
	auto xout = absge(c->posx[j] - posi[0], xmin);
	auto yout = absge(c->posy[j] - posi[1], xmin);
	auto zout = absge(c->posz[j] - posi[2], xmin);
	auto dx = c->cmposx[j] -posi[0];
	auto dy = c->cmposy[j] -posi[1];
	auto dz = c->cmposz[j] -posi[2];
	auto r2 = dx*dx + dy*dy+ dz*dz;
	// if (r2* theta2 > c->l[j]*c->l[j] && (xout ||yout||zout)){
	//     well_separated[j]=1;
	// } else{
	//     well_separated[j]=0;
	// }
	well_separated[j]=(r2* theta2 > c->l[j]*c->l[j] && (xout ||yout||zout));
    }
    for (auto j=0; j< NCHILDLEN; j++){
	if (well_separated[j]){
	    if (!c->isleaf[j] || c->nparticle[j]){
		*(pos_list+nlist) = vector(c->cmposx[j],
					   c->cmposy[j],
					   c->cmposz[j]);
		*(mass_list+nlist) = c->cmmass[j];
		nlist ++;
		if (nlist > list_max){
		    cerr << "List length exceeded\n";
		    exit(1);
		}
	    }
	}else{
	    // need to see if current node is this node.
	    bool self_interaction = false;
	    bhparticle * bp = bn + c->pindex[j];
	    //	    if ((b == bp)&&
	    //		li == c->l[j]){
	    if (c->posx[j]== posi[0]&&
		c->posy[j]== posi[1]&&
		c->posz[j]== posi[2] &&
		c->nparticle[j]){
		first_leaf = nlist;
		self_interaction= true;
	    }
	    if (self_interaction || c->isleaf[j]){
		for(auto i = 0; i < c->nparticle[j]; i++){
		    *(pos_list+nlist) = (bp+i)->get_rp()->get_pos();
		    
		    *(mass_list+nlist) =(bp+i)->get_rp()->get_mass();
		    nlist ++;
		    if (nlist > list_max){
			cerr << "List length exceeded\n";
			exit(1);
		    }
		}
	    }else{
		cn_add_to_interaction_list(cn,c->cindex[j],bn,
					   dest_node, theta2,
					   pos_list, mass_list,
					   nlist, list_max,
					   first_leaf);
	    }
	}
    }
}

void cn_bfs_add_to_interaction_list_step(child_nodes * cn,
					 int * index_list,
					 int * newindex,
					 int & nnodes,
					 bhparticle * bn,
					 bhnode & dest_node,
					 
					 real theta2,
					 vector * pos_list,
					 real * mass_list,
					 int & nlist,
					 int list_max,
					 int & first_leaf)
{
    //    printf("cn bfs step called with nnodes=%d\n", nnodes);
    int64_t well_separated[NCHILDLEN*nnodes];
    vector posi = dest_node.get_pos();
    real   li   = dest_node.get_length();
    bhparticle* b = dest_node.bpfirst;
    int newnnodes = 0;
    //#pragma fj loop unroll_count(8)    
    for(auto jn=0; jn<nnodes; jn++){
	auto c = cn+ index_list[jn];
	//	cn_opening_criterion(c,posi,li,theta2,well_separated+jn*NCHILDLEN);
	auto j8 = jn*NCHILDLEN;
#pragma omp simd
	//#pragma clang loop vectorize_width(8,scalable)
	//#pragma clang loop vectorize(enable)
	for (auto j=0; j< NCHILDLEN; j++){
	    auto jf = j8+j;
	    auto  xmin = (c->l[j] + li)*0.499999999999999;
	    auto xout = absge(c->posx[j] - posi[0], xmin);
	    auto yout = absge(c->posy[j] - posi[1], xmin);
	    auto zout = absge(c->posz[j] - posi[2], xmin);
	    auto dx = c->cmposx[j] -posi[0];
	    auto dy = c->cmposy[j] -posi[1];
	    auto dz = c->cmposz[j] -posi[2];
	    auto r2 = dx*dx + dy*dy+ dz*dz;
	    well_separated[jf]=(r2* theta2 > c->l[j]*c->l[j] && (xout ||yout||zout));
	}
    }
    
    for (auto j=0; j< NCHILDLEN*nnodes; j++){
	auto jn = j>>3;
	auto jc = j & (NCHILDLEN-1);
	auto c = cn+ index_list[jn];	
	if (well_separated[j]){
	    if (!c->isleaf[jc] || c->nparticle[jc]){
		*(pos_list+nlist) = vector(c->cmposx[jc],
					   c->cmposy[jc],
					   c->cmposz[jc]);
		*(mass_list+nlist) = c->cmmass[jc];
		nlist ++;
		if (nlist > list_max){
		    cerr << "List length exceeded\n";
		    exit(1);
		}
	    }
	}else{
	    // need to see if current node is this node.
	    bool self_interaction = false;
	    bhparticle * bp = bn + c->pindex[jc];
	    //	    if ((b == bp)&&
	    //		li == c->l[j]){
	    if (c->posx[jc]== posi[0]&&
		c->posy[jc]== posi[1]&&
		c->posz[jc]== posi[2] &&
		c->nparticle[jc]){
		first_leaf = nlist;
		self_interaction= true;
	    }
	    if (self_interaction || c->isleaf[jc]){
		for(auto i = 0; i < c->nparticle[jc]; i++){
		    *(pos_list+nlist+i) = (bp+i)->get_rp()->get_pos();
		    *(mass_list+nlist+i) =(bp+i)->get_rp()->get_mass();

		}
		nlist +=c->nparticle[jc] ;
		if (nlist > list_max){
		    cerr << "List length exceeded\n";
		    exit(1);
		}
	    }else{
		newindex[newnnodes]=c->cindex[jc];
		newnnodes ++;
	    }
	}
    }
    nnodes=newnnodes;
    for(int i=0;i<newnnodes; i++)index_list[i] = newindex[i];

}

void cn_bfs_add_to_interaction_list(child_nodes * cn,
				bhparticle * bn,
				bhnode & dest_node,

				real theta2,
				vector * pos_list,
				real * mass_list,
				int & nlist,
				int list_max,
				int & first_leaf)
{
    int nnodes=1;
    int index_list[list_max];
    int newindex[list_max];
    index_list[0] = 0;
    while (nnodes){
	cn_bfs_add_to_interaction_list_step(cn, index_list,  newindex,  nnodes,
					    bn, dest_node, theta2,
					    pos_list, mass_list, nlist,
					    list_max,first_leaf);
    }
}


void bhnode::bfs_add_to_interaction_list_step(bhnode * nodes,  real theta2,
					      int * index_list,
					      int * newindex,
					      int & nnodes,
					      vector * pos_list,
					      real * mass_list,
					      int & nlist,
					      int list_max,
					      int & first_leaf)
{
    int next_n=0;
    //    fprintf(stderr, "step called with nnodes= %d\n", nnodes);
    for(int i=0; i< nnodes; i++){
	//	fprintf(stderr, "i= %d\n", i);
	if (index_list[i] <0) continue;
	auto b = nodes+index_list[i];
	//	fprintf(stderr, "index= %d\n", index_list[i]);
	//	b->dump();
	auto l2 = b->l*b->l;
	if((separation_squared(this,b->cmpos)*theta2 > l2)
	   && (!are_overlapped(b,this) ) ){
	    // node and position is well separated;
	    //	    fprintf(stderr, "node %d added to list\n", index_list[i]);
	    *(pos_list+nlist) = b->cmpos;
	    *(mass_list+nlist) = b->cmmass;
	    nlist ++;
	    if (nlist > list_max){
		cerr << "List length exceeded\n";
		exit(1);
	    }
	}else{
	    int i;
	    if (b->isleaf || (b == (this))){
		if (b == (this)){
		    // adding the particles in the node itself
		    first_leaf = nlist;
		}
		bhparticle * bp = b->bpfirst;
		for(auto ii = 0; ii < b->nparticle; ii++){
 //		    fprintf(stderr, "ii = %d  nlist= %d\n", ii, nlist);
		    *(pos_list+nlist) = (bp+ii)->get_rp()->get_pos();
		    
		    *(mass_list+nlist) =(bp+ii)->get_rp()->get_mass();
		    nlist ++;
		    if (nlist > list_max){
			cerr << "List length exceeded\n";
			exit(1);
		    }
		}
	    }else{
		for(i=0;i<8;i++){
		    if (b->child[i] != NULL){
			newindex[next_n] = b->child[i] -nodes;
		    }else{
			newindex[next_n] = -1;
		    }
			
			// fprintf(stderr, "i= %d %d  newindex=%d\n",
			// 	i, next_n, newindex[next_n]);
		    next_n ++;
		}
	    }
	}
    }
    nnodes = next_n;
    for(int i=0;i<next_n; i++)index_list[i] = newindex[i];
}

void bhnode::bfs_add_to_interaction_list(bhnode * nodes,  real theta2,
					 bhnode * top,
					 vector * pos_list,
					 real * mass_list,
					 int & nlist,
					 int list_max,
					 int & first_leaf)
{
    int nnodes=1;
    int index_list[list_max];
    int newindex[list_max];
    index_list[0] = top - nodes;
    //    fprintf(stderr,"enter BFS_add_to_interaction_list\n");
    //    dump();
    while (nnodes){
	bfs_add_to_interaction_list_step(nodes, theta2, index_list, newindex,
					 nnodes,
					 pos_list, mass_list, nlist,list_max,
					 first_leaf);
    }
}


void bhnode::add_to_neighbour_list(bhnode & dest_node, real cutoff,
				     vector * pos_list,
				     real * mass_list,
				     int & nlist,
				     int list_max,
				     int & first_leaf)
{

    //    PRL(cutoff);
    if(!are_overlapped_with_cutoff(this,&dest_node, cutoff)){
	    // node and position is well separated
	    // just skip this node
	
    }else{
	int i;
	bhparticle * bp = bpfirst;
	//	PRC(nparticle); PRL(nplimit);
	if (isleaf||nparticle < nplimit || (this == (&dest_node))){
	    //	if (isleaf || (this == (&dest_node))){
	    if (this == (&dest_node)){
		// adding the particles in the node itself
		first_leaf = nlist;
		for(i = 0; i < nparticle; i++){
		    *(pos_list+nlist) = (bp+i)->get_rp()->get_pos();
		    
		    *(mass_list+nlist) =(bp+i)->get_rp()->get_mass();
		    nlist ++;
		    if (nlist > list_max){
			cerr << "List length exceeded\n";
			exit(1);
		    }
		}
	    }else{
		// for(i = 0; i < nparticle; i++){
		//     // need to judge distance here
		//     if (are_overlapped_with_cutoff(&dest_node, bp+i, cutoff)){
		// 	*(pos_list+nlist) = (bp+i)->get_rp()->get_pos();
			
		// 	*(mass_list+nlist) =(bp+i)->get_rp()->get_mass();
		// 	nlist ++;
		// 	if (nlist > list_max){
		// 	    cerr << "List length exceeded\n";
		// 	    exit(1);
		// 	}
		//     }
		// }
		add_to_neighbour_list_leafcell(dest_node, cutoff,
					       pos_list, mass_list,
					       nlist, list_max,
					       first_leaf);
	    }
	}else{
	    for(i=0;i<8;i++){
		if (child[i] != NULL){
		    child[i]->add_to_neighbour_list(dest_node, cutoff,
						      pos_list, mass_list,
						      nlist, list_max,
						      first_leaf);
		}
	    }
	}
    }
}

void bhnode::add_to_neighbour_list_leafcell(bhnode & dest_node, real cutoff,
				     vector * pos_list,
				     real * mass_list,
				     int & nlist,
				     int list_max,
				     int & first_leaf)
{
    int i;
    //    cerr << "lefcell add called with nparticle=" << nparticle <<endl;
    bhparticle * bp = bpfirst;
    if (nlist +nparticle> list_max){
	cerr << "List length exceeded\n";
	exit(1);
    }
    for(i = 0; i < nparticle; i++){
	// need to judge distance here
	if (are_overlapped_with_cutoff(&dest_node, bp+i, cutoff)){
	    *(pos_list+nlist) = (bp+i)->get_rp()->get_pos();
	    
	    *(mass_list+nlist) =(bp+i)->get_rp()->get_mass();
	    nlist ++;
	}
    }
}

void bhnode::add_to_neighbour_list_leafcell2(bhnode & dest_node,
					     real cutoff,
					     vector * pos_list,
					     real * mass_list,
					     int & nlist,
					     int list_max,
					     int & first_leaf)
{
    int i;
    //    cerr << "lefcell add called with nparticle=" << nparticle <<endl;
    bhparticle * bp = bpfirst;
    if (nlist +nparticle> list_max){
	cerr << "List length exceeded\n";
	exit(1);
    }
    real xmin = dest_node.get_length()*0.4999999999999999+cutoff;
    real xmin2 =dest_node.get_length()*0.4999999999999999;
    vector dpos = dest_node.get_pos();
#define LEAFMAX 256    
    real xa[LEAFMAX];
    real ya[LEAFMAX];
    real za[LEAFMAX];
    int  flags[LEAFMAX];
#pragma omp simd    
    for(i = 0; i < nparticle; i++){
	vector ppos = (bp+i)->get_rp()->get_pos();
	xa[i] = ppos[0];
	ya[i] = ppos[1];
	za[i] = ppos[2];
	flags[i]=1;
    }
#if 1    
#pragma omp simd    
    for(i = 0; i < nparticle; i++){
	real dx =xa[i]-dpos[0];
	real dy =ya[i]-dpos[1];
	real dz =za[i]-dpos[2];
	real r2sum = 0;
	real d = fabs(dx) - xmin2;
	if (d > 0) r2sum += d*d;
	d = fabs(dy) - xmin2;
	if (d > 0) r2sum += d*d;
	d = fabs(dz) - xmin2;
	if (d > 0) r2sum += d*d;
	if (r2sum <	cutoff*cutoff) {
	    flags[i]=1;
	}
    }
    for(i = 0; i < nparticle; i++){
	if (flags[i]){
	    
	    *(pos_list+nlist) = (bp+i)->get_rp()->get_pos();
	    
	    *(mass_list+nlist) =(bp+i)->get_rp()->get_mass();
	    nlist ++;
	}
    }
#else    
    for(i = 0; i < nparticle; i++){
	*(pos_list+nlist) = (bp+i)->get_rp()->get_pos();
	
	*(mass_list+nlist) =(bp+i)->get_rp()->get_mass();
	    nlist ++;
    }
#endif    
    //    nlist += nparticle;
}




#ifdef HARP3

extern "C" void h3open_();
extern "C" void h3close_();
extern "C" void accel_by_harp3_separate_noopen_(int * ni, vector * xi,
						int * nj, vector * xj,
						real *m,
						vector *  a,
						real *p,
						real * eps2);

void calculate_force_from_interaction_list_using_grape4(vector * pos_list, real * mass_list,
							int list_length, int first_leaf, int ni,
							real eps2,
							vector * acc_list, real * phi_list)
{
    static call_count = 0;
    static h3_open_state = 0;
    if (h3_open_state == 0){
	h3open_();
	h3_open_state = 1;
    }
    //    PR(ni);PRL(list_length);
    nisum += ni;
    tree_walks += 1;
    total_interactions += ((real)ni)*list_length;
    accel_by_harp3_separate_noopen_(&ni,pos_list+first_leaf, &list_length,pos_list, mass_list,
				    acc_list, phi_list, &eps2);
    call_count += ni;
    if (call_count > 500000){
	cerr << "Close and release GRAPE-4\n";
	h3close_();
	h3_open_state = 0;
	call_count = 0;
    }
}
#endif

void bhnode::evaluate_gravity_using_tree_and_list(bhnode & source_node,
						  child_nodes * cn,
						  bhparticle * bn,
						  real theta2,
						  real eps2,
						  int ncrit)
{
    const int list_max = LISTMAX;
    static real mass_list[list_max];
    static vector pos_list[list_max];
    real epsinv = 1.0/sqrt(eps2);

    static vector * acc_list = NULL;
    static real * phi_list = NULL;
    if (acc_list == NULL){
	acc_list = new vector[ncrit + 100];
	phi_list = new real[ncrit + 100];
    }

    //    PR(pos); PR(nparticle); PRL(isleaf);
    if((nparticle > ncrit) && (isleaf==0)){
	for(int i=0;i<8;i++){
	    if (child[i] != NULL){
		child[i]->evaluate_gravity_using_tree_and_list(source_node,
							       cn, bn,
							       theta2,
							       eps2,
							       ncrit);
	    }
	}
    }else{
	//
	// node is below critical ... first create list
	//
	int list_length = 0;
	int first_leaf = -1;
	double toff = NbodyLib::GetWtime();
#if defined(DFS_TREE_WALK)
	source_node.add_to_interaction_list(*this, theta2,
					    pos_list,
					    mass_list,
					    list_length,
					    list_max,
					    first_leaf);
#elif defined(SIMD_DFS_TREE_WALK)
	cn_add_to_interaction_list(cn, 0, bn,
				   *this, theta2,
				   pos_list,
				   mass_list,
				   list_length,
				   list_max,
				   first_leaf);
#elif defined(SIMD_TREE_WALK)
	cn_bfs_add_to_interaction_list(cn,  bn,
				       *this, theta2,
				       pos_list,
				       mass_list,
				       list_length,
				       list_max,
				       first_leaf);
#else
	bfs_add_to_interaction_list(&source_node,  theta2, &source_node,
				    pos_list,  mass_list,
				    list_length,
				    list_max,
				    first_leaf);
#endif	
	if (first_leaf == -1){
	    cerr << "evaluate_gravity: impossible error \n";
	    cerr << "failed toC find the node in the tree \n";
	    exit(1);
	}
	
	bhparticle * bp = bpfirst;

#if 0
	printf("nparticle= %d nlist= %d fl= %d\n", nparticle, list_length, first_leaf);
	for(auto j=0;j<list_length; j++){
	    printf("%d: %g %g %g %g\n",j, pos_list[j][0],
		   pos_list[j][1],pos_list[j][2], mass_list[j]);
	}
#endif	
#ifndef HARP3
	tree_walks ++;
	nisum += nparticle;
	total_interactions += ((real)nparticle)*list_length;
	double toff2=NbodyLib::GetWtime();
	t_walk += toff2  - toff;
#pragma omp parallel for
	for(int i = 0; i < nparticle; i++){
	    real_particle * p = (bp+i)->get_rp();
	    vector acc;
	    real phi;
	    calculate_force_from_interaction_list(pos_list[i+first_leaf],eps2, acc, phi,
					  pos_list,mass_list,list_length);
	    p->set_acc_gravity(acc);
	    p->set_phi_gravity(phi + p->get_mass()*epsinv);
	}
	double toff3=NbodyLib::GetWtime();
	t_calc += toff3  - toff2;
#else
	calculate_force_from_interaction_list_using_grape4(pos_list, mass_list,list_length, first_leaf,
							   nparticle, eps2, acc_list, phi_list);
	for(int i = 0; i < nparticle; i++){
	    real_particle * p = (bp+i)->get_rp();
	    p->set_acc_gravity(acc_list[i]);
	    p->set_phi_gravity(phi_list[i] + p->get_mass()*epsinv);
	}
#endif	
    }
}


void bhnode::make_neighbour_list_using_tree(bhnode & source_node,
					    real cutoff,
					    int ncrit)
{
    const int list_max = 40000;
    static real mass_list[list_max];
    static vector pos_list[list_max];
    static vector * acc_list = NULL;
    static real * phi_list = NULL;
    if (acc_list == NULL){
	acc_list = new vector[ncrit + 100];
	phi_list = new real[ncrit + 100];
    }
    //    PR(pos); PR(nparticle); PRL(isleaf);
    if((nparticle > ncrit) && (isleaf==0)){
	for(int i=0;i<8;i++){
	    if (child[i] != NULL){
		child[i]->make_neighbour_list_using_tree(source_node,
							 cutoff,
							 ncrit);
	    }
	}
    }else{
	//
	// node is below critical ... first create list
	//
	int list_length = 0;
	pilist = new interaction_list[1];
	pilist->ni = nparticle;
	source_node.add_to_neighbour_list(*this,  cutoff,
					    pilist->pos_list,
					    pilist->mass_list,
					    pilist->length,
					    list_max,
					    pilist->first_leaf);
	if (pilist->first_leaf == -1){
	    cerr << "evaluate_gravity: impossible error \n";
	    cerr << "failed to find the node in the tree \n";
	    exit(1);
	}
	list_length=pilist->length;
	bhparticle * bp = bpfirst;
	tree_walks ++;
	nisum += nparticle;
	total_interactions += ((real)nparticle)*list_length;
    }
}

void bhnode::make_neighbour_list_using_dual_treewalk(bhnode & source_node,
						     real cutoff,
						     int ncrit,
						     int srclevel,
						     int destlevel)
{
    //    PRC(srclevel),PRL(destlevel);
    //    PR(pos); PR(nparticle); PRL(isleaf);
    int src_down = 0;
    int dest_down = 0;
    if(!are_overlapped_with_cutoff(this,&source_node, cutoff)){
	// node and position is well separated
	// just skip this node
	return;
    }
    if(nparticle > ncrit){
	if (source_node.is_leaf()||source_node.nparticle < nplimit){
	    dest_down = 1;
	}else{
	    if (srclevel < destlevel){
		src_down = 1;
	    }else{
		dest_down = 1;
	    }
	}
    }else{
	if (pilist == NULL) {
	    pilist = new interaction_list[1];
	    pilist->ni = nparticle;
	}
	//	if ((!source_node.is_leaf())
	if ((!source_node.is_leaf() &&(source_node.nparticle < nplimit ))
	    &&(this != &source_node)) src_down = 1;
    }
    if(src_down){
	for(int i=0;i<8;i++){
	    if (source_node.child[i] != NULL){
		make_neighbour_list_using_dual_treewalk(*(source_node.child[i]),
							cutoff,
							ncrit,
							srclevel+1,
							destlevel);
	    }
	}
    }else if (dest_down){
	for(int i=0;i<8;i++){
	    if (child[i] != NULL){
		child[i]->make_neighbour_list_using_dual_treewalk(source_node,
								  cutoff,
								  ncrit,
								  srclevel,
								  destlevel+1);
	    }
	}
	
    }else{
	//
	// node is below critical ... first create list
	//
	source_node.add_to_neighbour_list(*this,  cutoff,
					    pilist->pos_list,
					    pilist->mass_list,
					    pilist->length,
					    ilist_max,
					    pilist->first_leaf);
    }
}


void bhnode::dump_interaction_list()
{
    if(pilist != NULL){
	pilist->dump();
    }else{
	for(int i=0;i<8;i++){
	    if (child[i] != NULL){
		child[i]->dump_interaction_list();
	    }
	}
    }
}
void bhnode::clear_interaction_list()
{
    if(pilist != NULL){
	//	delete [] pilist;
	//	pilist = NULL;
	pilist-> length =0;
    }else{
	for(int i=0;i<8;i++){
	    if (child[i] != NULL){
		child[i]->clear_interaction_list();
	    }
	}
    }
}
void bhnode::update_tree_counters()
{
    if(pilist != NULL){
	total_interactions += ((real)(pilist->ni))*pilist->length;
	tree_walks ++;
	nisum += pilist->ni;
    }else{
	for(int i=0;i<8;i++){
	    if (child[i] != NULL){
		child[i]->update_tree_counters();
	    }
	}
    }
}

void real_system::evaluate_gravity_using_default_tree_and_list(real theta2,
					  real eps2,
					  int ncrit)
{
    //    bn->dump();
    bn->evaluate_gravity_using_tree_and_list(*bn, cn, bp,  theta2,eps2, ncrit);
}

void real_system::make_neighbour_list_using_default_tree(real cutoff,
					    int ncrit,int nctree)
{
    cerr << "Enter make neighbour list , cpu = " <<NbodyLib::GetWtime() << endl;
    double tstart = NbodyLib::GetWtime();
    bn->set_nplimit(nctree);
    bn->make_neighbour_list_using_tree(*bn, cutoff, ncrit);
    cerr << "Exit make neighbour list , dt = " <<NbodyLib::GetWtime()-tstart << endl;
    //    bn->dump_interaction_list();
}
void real_system::make_neighbour_list_using_default_tree_and_dual_walk(real cutoff,
							  int ncrit,
							  int nctree)
{
    clear_tree_counters();
    bn->set_nplimit(nctree);
    bn->clear_interaction_list();
    cerr << "Enter make neighbour list DT, cpu = " <<NbodyLib::GetWtime() << endl;
    double tstart = NbodyLib::GetWtime();
    bn->make_neighbour_list_using_dual_treewalk(*bn, cutoff, ncrit,0,0);
    cerr << "Exit make neighbour list DT, dt = " <<NbodyLib::GetWtime()-tstart << endl;
    //    bn->dump_interaction_list();
    bn->update_tree_counters();
}

// void real_particle::calculate_gravity_using_tree(real eps2, real theta2)
// {
//     acc_gravity = 0;
//     phi_gravity = mass/sqrt(eps2);
//     bn->accumulate_force_from_tree(pos,eps2,theta2,
// 				  acc_gravity, phi_gravity);
//     nisum += 1;
// }


#ifdef TESTXXX
//
// do the serious test of
// construction of tree
// consistency of tree
// validity of the neighbour list (by comparing with the result
// of direct calculation

void main()
{
    static sph_system pb;
    int n;
    cerr << "Enter n:";
    cin >> n ;
    pb.create_uniform_sphere(n, 0 , 1);
    //    pb.dump();
    int nkey = 0;
    pb.initialize_h_and_nbl(pow(1.0/n,0.33333));
    static sph_system pbcopy = pb;
    copy_sph_particles(&pb, &pbcopy);
    real_particle * psph = pb.get_particle_pointer();
    real_particle * psphcopy = pbcopy.get_particle_pointer();
    pb.set_nnb_using_tree();
    cerr << "Dumping copy ... \n";
    cerr << "checking NB \n";
    int error = 0;
    for(int i = 0; i<n; i++){
	(psph+i)->sort_nblist();
	int err = 0;
	if((psph+i)->get_nnb() != (psphcopy+i)->get_nnb()){
	    cerr << "Neighbour count differs for "; PRL(i);
	    err = 1;
	    
	}
	if (err == 0){
	    for(int j = 0; (j< (psph+i)->get_nnb()) && (err == 0); j++){
		if ((psph+i)->get_neighbour(j)->get_index()!=
		    (psphcopy+i)->get_neighbour(j)->get_index()) err = 1;
	    }
	}
	if(err){
	    (psph+i)->dump();
	    (psphcopy+i)->dump();
	    error ++;
	}
    }
    PRL(error);
}
#endif

#ifdef TEST
//
// do the serious test of
// construction of tree
// consistency of tree
// validity of the neighbour list (by comparing with the result
// of direct calculation

void main()
{
    static sph_system pb;
    int n;
    cerr << "Enter n:";
    cin >> n ;
    pb.create_uniform_sphere(n, 0 , 1);
    //    pb.dump();
    int nkey = 0;
    bhparticle * bp = NULL;

    bn = new bhnode[n];
    for(int i = 0; i<1; i++){
	real rsize = initialize_key(n,pb.get_particle_pointer(),nkey,bp);
	for(int j = 0; j<n;j++) (bn+j)->clear();
	bn->assign_root(vector(0.0), rsize*2, bp, n);
	bhnode * btmp = bn+1;
	int heap_remainder = n-1;
	BHlong key = 0;
        bn->create_tree_recursive(btmp,  heap_remainder,key, default_key_length, 4);
    }
    PRL(bn->sanity_check());
    pb.initialize_h_and_nbl(pow(1.0/n,0.33333));
    bn->set_hmax_for_sph();
    //    bn->dump();
    bn->set_cm_quantities();
    //    bn->dump();
    static sph_system pbcopy = pb;
    copy_sph_particles(&pb, &pbcopy);
    real_particle * psph = pb.get_particle_pointer();
    real_particle * psphcopy = pbcopy.get_particle_pointer();
    for(int i = 0; i<n; i++){
	(psph+i)->clear_nnb();
    }
    PRL(check_and_set_nbl(bn, bn));
    cerr << "Dumping copy ... \n";
    cerr << "checking NB \n";
    int error = 0;
    for(int i = 0; i<n; i++){
	(psph+i)->sort_nblist();
	int err = 0;
	if((psph+i)->get_nnb() != (psphcopy+i)->get_nnb()){
	    cerr << "Neighbour count differs for "; PRL(i);
	    err = 1;
	    
	}
	if (err == 0){
	    for(int j = 0; (j< (psph+i)->get_nnb()) && (err == 0); j++){
		if ((psph+i)->get_neighbour(j)->get_index()!=
		    (psphcopy+i)->get_neighbour(j)->get_index()) err = 1;
	    }
	}
	if(err){
	    (psph+i)->dump();
	    (psphcopy+i)->dump();
	    error ++;
	}
    }
    PRL(error);
    pb.use_self_gravity = 1;
    pb.eps2_for_gravity = 0.01;
#define COMPARISON_WITH_DIRECT    
#ifdef COMPARISON_WITH_DIRECT    
    pb.calculate_uncorrected_gravity_direct();
    copy_sph_particles(&pb, &pbcopy);
    psphcopy = pbcopy.get_particle_pointer();
    cerr << "Direct force \n";
    for(int i = 0; i<n; i++){
	real phi = (psphcopy+i)->get_phi_gravity();
	vector acc  = (psphcopy+i)->get_acc_gravity();
	PR(i); PR(phi); PRL(acc);
    }
#endif


    
    cerr << "Tree   force \n";
    for(int j = 0; j<10; j++){
	PRL(j);
	pb.apply_vf(real_particle::clear_acc_phi_gravity);
	for(int i = 0; i<n; i++){
	    (psph+i)->calculate_gravity_using_tree(pb.eps2_for_gravity, 0.4);
	}
    }
    pb.apply_vf(real_particle::clear_acc_phi_gravity);
    bn->evaluate_gravity_using_tree_and_list(*bn,0.4,pb.eps2_for_gravity,1);
#ifdef COMPARISON_WITH_DIRECT    
    real perrmax = 0;
    real ferrmax = 0;
    for(int i = 0; i<n; i++){
	real phi = (psph+i)->get_phi_gravity();
	real phierr = (psphcopy+i)->get_phi_gravity()-phi;
	vector acc  = (psph+i)->get_acc_gravity();
	vector accerr  = (psphcopy+i)->get_acc_gravity()-acc;
	PR(i); PR(phi); PRC(acc); PRC(phierr); PRL(accerr);
	real prelerr = fabs(phierr/phi);
	real frelerr = abs(accerr)/abs(acc);
	if(perrmax < prelerr) perrmax = prelerr;
	if(ferrmax < frelerr) ferrmax = frelerr;
    }
    PR(perrmax);    PRL(ferrmax);
#else
    for(int i = 0; i<n; i++){
	real phi = (psph+i)->get_phi_gravity();
	vector acc  = (psph+i)->get_acc_gravity();
	PR(i); PR(phi); PRL(acc); 
    }
	
#endif    
    
}
#endif

	

#ifdef TESTXX
//
// Sample test for timing purpose...

void main()
{
    static sph_system pb;
    int n;
    cerr << "Enter n:";
    cin >> n ;
    pb.create_uniform_sphere(n, 0 , 1);
    //    pb.dump();
    int nkey = 0;
    bhparticle * bp = NULL;

    bn = new bhnode[n];
    for(int i = 0; i<10; i++){
	real rsize = initialize_key(n,pb.get_particle_pointer(),nkey,bp);
	for(int j = 0; j<n;j++) (bn+j)->clear();
	bn->assign_root(vector(0.0), rsize*2, bp, n);
	bhnode * btmp = bn+1;
	int heap_remainder = n-1;
	BHlong key = 0;
        bn->create_tree_recursive(btmp,
				  heap_remainder,key,
				  default_key_length, 8 );
	PRL(heap_remainder);
    }
    PRL(bn->sanity_check());
    real_particle * psph = pb.get_particle_pointer();
    real h0 = pow(1.0/n,0.33333);
    for(int i = 0; i<10; i++){
	
	pb.apply_vf(real_particle::set_h, h0);
	pb.apply_vf(real_particle::clear_nnb);
	bn->set_hmax_for_sph();
	//    bn->dump();
	PRL(check_and_set_nbl(bn, bn));
	pb.apply_vf(real_particle::sort_nblist);
		
    }
}
#endif






