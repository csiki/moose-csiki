/**********************************************************************
** This program is part of 'MOOSE', the
** Messaging Object Oriented Simulation Environment.
**   copyright (C) 2003-2007 Upinder S. Bhalla, Niraj Dudani and NCBS
** It is made available under the terms of the
** GNU Lesser General Public License version 2.1
** See the file COPYING.LIB for the full notice.
**********************************************************************/

#ifndef _HSOLVE_PASSIVE_H
#define _HSOLVE_PASSIVE_H
#include "../basecode/header.h"
#include "../biophysics/CompartmentBase.h"
#include "../biophysics/Compartment.h"
using namespace moose; // For moose::Compartment from 'Compartment.h'
#include "HSolveUtils.h"
#include "HSolveStruct.h"
#include "../basecode/SparseMatrix.h"
#include "../diffusion/FastMatrixElim.h"
#include <map>
#include <set>
#include <algorithm> // std::set_intersection
#include <iterator> // std::inserter, std::advance

struct TreeNodeStruct
{
    set< unsigned int > children;	///< Hines indices of child compts
    // a set would be more
    double Ra;
    double Rm;
    double Cm;
    double Em;
    double initVm;
};

class HSolvePassive
{
#ifdef DO_UNIT_TESTS
	friend void testHSolvePassive();
#endif
	
public:
	void setup( Id seed, double dt );
	void solve();
	
protected:
	// Integration
	void updateMatrix();
	
    /**
     * Own stuff TODO
    */
    
    /**
     * Diagonal values of the Hines matrix. Needed to be saved separately
     * from the matrix, to initialise the matrix at every step.
     * Ordered according to Hines numbering.
     */
    vector< double > diagvals_;
    
    /**
     * Offdiagonal values of the Hines matrix. Needed to be saved separately
     * from the matrix, to initialise the matrix at every step.
     * The key of the map indicates the row and col of the matrix (Hines indices),
     * respectively. It's not redundant so there's no entry for (1,2)
     * and (2,1) connections. The value of the map is Gij = Gi * Gj / Gsum
     * taking into consideration junctions more than 2 compartments.
     */
    map< pair< unsigned int, unsigned int >, double > junctions_;
    
    /**
     * FastMxElim attributes TODO
    */
    FastMatrixElim passiveElim_;
    vector< Triplet<double> > passiveOps_;
    vector< double > passiveDiagVal_;
    
    /**
     * Attributes from HinesMatrix. TODO !!!
    */
    unsigned int              nCompt_;
    vector< double >          VMid_;    ///< Compartment voltage at the
    ///< middle of a time step.
    double                    dt_;
    
	vector< CompartmentStruct >       compartment_;
	vector< Id >                      compartmentId_;
	vector< double >                  V_;				/**< Compartment Vm.
		* V_ is addressed using a compartment index. V_ stores the Vm value
		* of each compartment. */
	vector< TreeNodeStruct >          tree_;			/**< Tree info.
		* The tree is used to acquire various values during setup. It contains
		* the user-defined original values of all compartment parameters.
		* Therefore, it is also used during reinit. */
	map< unsigned int, InjectStruct > inject_;			/**< inject map.
		* contains the list of compartments that have current injections into
		* them. */
	
private:
    
    // originally in HinesMx - TODO
    void prepareSparseMatrix();
    
	// Setting up of data structures
	void clear();
	void walkTree( Id seed );
	void initialize();
	void storeTree();
	
	// Used for unit tests.
	double getV( unsigned int row ) const;
};

#endif // _HSOLVE_PASSIVE_H
