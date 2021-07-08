/**
 *  txtfileio.h
 *
 *  A class to input/output data to text files.
 *
 *  Created by Yinan Li on August 09, 2018.
 *
 *  Hybrid Systems Group, University of Waterloo.
 */


#ifndef _txtfileio_h
#define _txtfileio_h

#include <fstream>
#include <iostream>
#include <vector>

#include "grid.h"
#include "csolver.h"
#include "transition.hpp"
#include "dsolver.h"
#include "buchi.h"


namespace rocs {

    /**
     * A txt file writer class.
     */
    class txtFileHandler {
    private:
	const char* _filename;
	std::fstream _txtfile;

    public:
	txtFileHandler(const char* f) : _filename(f) {}

	/**
	 * Open a txt file.
	 **/
	bool open(std::ios_base::openmode mode) {
	    _txtfile.open(_filename, mode);
	    return _txtfile.is_open();
	}
	/**
	 * Open a file by name.
	 * @param[in] fname the file name.
	 * @param[in] mode the output mode.
	 */
	bool open(const char* fname, std::ios_base::openmode mode) {
	    _txtfile.open(fname, mode);
	    return _txtfile.is_open();
	}
    
	/**
	 * Close a txt file.
	 **/
	void close() {_txtfile.close();}

	/**
	 * Save grid data to .txt file.
	 * @param filename the saved file name.
	 */
	void write_uniform_grids(const grid &g);

	/**
	 * Write leaves under a node to a .txt file (Depth-first).
	 * @param ctlr the controller in the form of SPtree.
	 * @param ptrn pointer to the given node.
	 */
	void write_sptree_leaves(const SPtree &ctlr, SPnode *ptrn);

	/**
	 * Write the controller generated by CSolver into a .txt file.
	 *
	 * The controller is saved in 3 variables:
	 * + pavings: the interval boxes,
	 * + ctlr: a list of binary values indicating the admissible control inputs, and
	 * + tags: the tags (0,1,or 2) of the boxes.
	 *
	 * @param sol the solver which contains the controller info.
	 **/
	void write_sptree_controller(const CSolver &sol) {
	    if (!_txtfile.is_open()) {open(std::ios::out);}
	    write_sptree_leaves(sol._ctlr, sol._ctlr._root);
	}

	/**
	 * Write the post transitions of an abstraction within a certain area into a txt file.
	 *
	 * Write all the transitions can be very time consuming.
	 * In the format of:
	 * #: (from state id) (control id) (to state id)
	 *
	 * @param nts the transition system.
	 * @param g the grid of states for the abstraction.
	 * @param xlb the lower bound of the given area.
	 * @param xub the upper bound of the given area.
	 **/
	void write_post_transitions(const fts &nts, const grid &g,
				    const double xlb[], const double xub[]);
	void write_pre_transitions(const fts &nts, const grid &g,
				    const double xlb[], const double xub[]);


	/**
	 * Read a discrete controller generated by a buchi game solver from a txt file.
	 *
	 * @param w_x0 an array of winning states.
	 * @param encode3 indexing array.
	 * @param nts_ctrlr an array of NODE_POST for looking up the ctrl pair.
	 * @param ctrl an array of (q,x) pair.
	 * @param q_prime an array storing the automaton transitions.
	 **/
	void read_discrete_controller(std::vector<long long> &w_x0,
				      std::vector<long long> &encode3,
				      std::vector<NODE_POST> &nts_ctrlr,
				      std::vector<CTRL> &ctrl,
				      std::vector<int> &q_prime);
    
    };//txtFileHandler class

}//namespace rocs


#endif
