#ifndef  BHTREE_H
#  define   BHTREE_H
/*-----------------------------------------------------------------------------
 *  BHtree : basic class for C++ implementation of BH treecode
 *  J. Makino 1998/12/14
 *-----------------------------------------------------------------------------
 */

    

const int default_key_length = 15;
const int default_ix_offset = 1<<(default_key_length-1);
// First include some particle class
// This need, someday, to be changed to template...
#include "vector.h"
#include "nbody_particle.h"
#define NBODY
#ifdef NBODY
typedef nbody_particle real_particle;
typedef nbody_system real_system;
#endif
#ifdef SPH
typedef sph_particle real_particle;
typedef sph_system real_system;
#endif

// The following dependence on ALPHA is introduced by
// Andreas.Adelmann@psi.ch
// Thank you Andreas!

#ifdef __alpha
   typedef  long BHlong;
#else
   typedef  long long BHlong;
#endif


const int ilist_max = 1000;

class  interaction_list{
 public:
    int length;
    int first_leaf;
    int ni;
    real mass_list[ilist_max];
    vector pos_list[ilist_max];
    interaction_list(){
	length=0;
	first_leaf=-1;
	ni = 0;
    }
    void dump()
    {
	cout <<  "Neighbours " << length <<endl;
	for(int i = 0; i < length; i++){
	    cout <<  pos_list[i] <<endl;
	}
	cout << "i particles " << ni << endl;
	for (int i=0;i<ni;i++){
	    cout << pos_list[first_leaf+i]<< endl;
	}
    }

	    
};

class bhparticle
{
private:
    real_particle * rp;
    BHlong key;
    
public:
    bhparticle(){
	rp = NULL;
	key = 0;
    }
    void set_key(real rscale, int ioffset, int keybits);
    int friend compare_key( bhparticle * p1,  bhparticle * p2);
    BHlong get_key(){return key;}
    void set_rp(real_particle * p){rp = p;}
    real_particle * get_rp(){return rp ;}
    void friend sort_bh_array( bhparticle * r, int lo, int up );
    
};


class bhnode
{
private:
    static int nplimit;
    vector pos;
    real l;
    bhnode * child[8];
    bhparticle * bpfirst;
    int nparticle;
    int isleaf;
#ifdef SPH    
    real hmax_for_sph;
#endif    
    vector cmpos;
    real cmmass;
    
public:
    interaction_list * pilist;
    bhnode(){
	pos = 0.0;
	l = 0.0;
	for(int i = 0; i<8;i++)child[i] = NULL;
	bpfirst = NULL;
	pilist=NULL;
	nparticle = 0;
	isleaf = 1;
#ifdef SPH	
	hmax_for_sph = 0;
#endif	
	cmpos = 0.0;
	cmmass = 0.0;
    }
    void clear(){
	pos = 0.0;
	l = 0.0;
	for(int i = 0; i<8;i++)child[i] = NULL;
	bpfirst = NULL;
	nparticle = 0;
	isleaf = 1;
#ifdef SPH	
	hmax_for_sph = 0;
#endif	
	cmpos = 0.0;
	cmmass = 0.0;
    }
    int is_leaf(){return isleaf;}
    void set_nplimit(int n){nplimit = n;}
    void set_pos(vector newpos){pos = newpos;}
    vector get_pos(){return pos;}
	vector*  get_posp()                         {return &pos;}
    void set_length(real newl){l = newl;}
    real get_length(){return l;}
    void create_tree_recursive(bhnode * & heap_top, int & heap_remainder,
			       BHlong current_key,
			       int current_level,
			       int n_critical);
    void assign_root(vector root_pos, real length, bhparticle * bp, int nparticle);
    void dump(int indent = 0);
    int sanity_check();
#ifdef SPH    
    void set_hmax_for_sph();
    real get_hmax_for_sph(){return hmax_for_sph;}
#endif    
    int friend check_and_set_nbl(bhnode * p1,bhnode * p2);
    void set_cm_quantities();
    void accumulate_force_from_tree(vector & ipos, real eps2, real theta2,
				   vector & acc,
				   real & phi);
    void add_to_interaction_list(bhnode & dest_node, real theta2,
				 vector * pos_list,
				 real * mass_list,
				 int & nlist,
				 int list_max,
				 int & first_leaf);
    void add_to_neighbour_list(bhnode & dest_node, real cutoff,
				      vector * pos_list,
				      real * mass_list,
				      int & nlist,
				      int list_max,
				      int & first_leaf);
    void evaluate_gravity_using_tree_and_list(bhnode & source_node,
					      real theta2,
					      real eps2,
					      int ncrit);
    void make_neighbour_list_using_tree(bhnode & source_node,
					real cutoff,
					int ncrit);

    void make_neighbour_list_using_dual_treewalk(bhnode & source_node,
						 real cutoff,
						 int ncrit,
						 int srclevel,
						 int destlevel);
    void dump_interaction_list();
    void clear_interaction_list();
    void update_tree_counters();
};


void clear_tree_counters();
void print_tree_counters();

    
#endif
