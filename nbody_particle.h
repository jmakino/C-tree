#pragma once
/*-----------------------------------------------------------------------------
 *  nbody-particle : basic class for simple nbody implementation
 *  J. Makino 1998/11/29
 *-----------------------------------------------------------------------------
 */
#include "vector.h"

#ifndef ONED
#define THREED
#endif
#define REAL_GRAVITY


class nbody_particle
    {


    public:
        vector pos;
	vector vel;
	vector acc_gravity;
	real phi_gravity;
	real phi_gravity_external;
	real mass;
	int index;
	nbody_particle(){
	    pos = 0.0;
	    vel = 0.0;
	    acc_gravity = 0.0;
	    phi_gravity =  mass =  0.0;
	    index = 0;
	}
        void  set_pos(const vector& new_pos)      {pos = new_pos;}
        void  set_vel(const vector& new_vel)      {vel = new_vel;}
        void  set_acc_gravity(const vector& new_acc)      {acc_gravity = new_acc;}
        void  set_phi_gravity(real new_phi)      {phi_gravity = new_phi;}

	void  clear_pos()                         {pos = 0.0;}
	void  clear_vel()                         {vel = 0.0;}
	void  clear_acc_phi_gravity(){acc_gravity = 0.0;phi_gravity = 0.0;}
	void  correct_phi_self_gravity(real epsinv)      {phi_gravity += mass*epsinv;}	

	void  inc_pos(const vector& d_pos)        {pos += d_pos; }
	void  inc_vel(const vector& d_vel)        {vel += d_vel;}
	void  update_vel(real dt)        {vel += dt*acc_gravity;}
	void  update_pos(real dt)        {pos = (pos+dt*vel).readjust();}

	void  scale_pos(const real scale_factor)  {pos *= scale_factor; }
	void  scale_vel(const real scale_factor)  {vel *= scale_factor; }

	vector  get_pos()                         {return pos;}
	vector*  get_posp()                         {return &pos;}

	vector  get_vel()                         {return vel;}
	real  get_phi_gravity()                         {return phi_gravity;}
	real  get_phi_gravity_external()  {return phi_gravity_external;}
	vector  get_acc_gravity()                         {return acc_gravity;}

	real get_mass()			{return mass;}
	void set_mass(real m)		{mass = m;}
	void set_index(int i){index = i;}
	int get_index(){return index;}
	void predict(real dt){
	    real dt2 = dt*dt*0.5;
	    pos = (pos + dt*vel + dt2*acc_gravity).readjust();
	    vel +=  (dt*0.5)*acc_gravity;
	}
	void correct(real dt){
	    vel +=  (dt*0.5)*acc_gravity;
	}

	void read(istream & );
	void write(ostream & );
	void dump();

	void plot(real parm);

	real kinetic_energy();
	real energy();
	real get_ke(){ return 0.5*mass*vel*vel;}

	void friend accumulate_mutual_gravity(nbody_particle & p1,
					      nbody_particle & p2,
					      real eps2);
	void calculate_gravity_using_tree(real eps2, real theta2);
	void correct_gravity(real);
	void external_field(real);

};

typedef vector (nbody_particle::*nbody_VMF_ptr)(void);     // vector member function pointer
typedef void (nbody_particle::*nbody_MF_ptr)(const vector &);     // member function pointer

typedef void (nbody_particle::*nbody_VF_ptr)(void);     // void member function pointer
typedef void (nbody_particle::*nbody_RF_ptr)(real);     // void member function
						    // pointer with real arg
typedef void (nbody_particle::*nbody_RRF_ptr)(real,real);     // void member function
						    // pointer with two real args
