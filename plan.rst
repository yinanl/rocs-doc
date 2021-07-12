Work In Progress
================

.. note::
   This section is kept as a project description of the work under development.


LTL-to-BA translator embedding
------------------------------
The fact that the user has to manually transcribe the DBA converted from an LTL specification into a text file in a specific form is not convenient. To improve the usability, we are working on integrating the LTL-to-BA translation function in Spot to ROCS.

How is a specification file parsed in ROCS
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
As we have shown in the :ref:`User Guide <sec_spec_file>`, the specification file includes two groups of information, and they are parsed by ``DBAparser``:

- The metadata that includes the number of DBA states (or nodes) ``nNodes``, the number of atomic propositions ``nAP``, the accepting nodes ``acc``, and the initial state ``q0``. These four pieces of information will be read directly, by a ``DBAparser``, to four local variables (defined by the user) in the main C++ program of each case study. These variables will be passed to the solvers for control synthesis.

- All the transitions in the DBA in the matrix form ``M``. Currently, a 2d array variable ``arrayM`` (of type ``std::vector<std::vector<rocs::UintSmall>>``) needs to be defined in the main program for temporarily holding the transitions parsed by the ``DBAparser`` before they are passed to the solvers (either ``CSolver`` or ``BSolver``). The following is a diagram to illustrate such a workflow.

.. graphviz::

   digraph parse {
       rankdir="LR";
       M [shape=box];
       M->arrayM [label="DBAparser"];
       arrayM [shape=box];
	   arrayM->_M [label="CSolver"];
	   _M [shape=box];
	   arrayM -> qprime [label="BSolver"];
	   qprime [shape=box, label="_sol._dba->q_prime"];
   }

``_M`` is a member variable of the ``CSolver`` class, which performs specification-guided synthesis, and ``sol`` is a member variable of the ``BSolver``, which is responsible for abstraction-based control synthesis. Since abstraction-based engine is coded in pure C, the ``_sol`` is a head that covers all the other data structures used in the C code. The ``_dba`` is a data structure that stores the DBA information in ``_sol``, which holds the pointer to the ``q_prime`` (that stores all the transitions of a DBA).

Each row of the 2d variable ``arrayM`` represents the ending states in the DBA under different propositions from each DBA node. This means that for a DBA state :math:`q` (the index of type ``rocs::UintSmall``), the DBA state to which the transition with the proposition :math:`p` (the decimal value of the binary encoding) leads is ``arrayM[q][p]``, which is also the index of a DBA state. The ``arrayM`` can be visualized as follows.


+---------------------+-------------+-------------+-----+
| node(r),prop(c) ids | 0           | 1           | ... |
+=====================+=============+=============+=====+
| 0                   | array[0][0] | array[0][1] | ... |
+---------------------+-------------+-------------+-----+
| 1                   | array[1][0] | ...         | ... |
+---------------------+-------------+-------------+-----+
| ...                 | ...         | ...         | ... |
+---------------------+-------------+-------------+-----+



Project decomposition
^^^^^^^^^^^^^^^^^^^^^
Embedding the LTL-BA translation by using Spot may consist:

- Simplify the DBA translated by Spot. This at least includes

  * trimming off the transitions that are impossible to be triggered, and
  * trimming off the isolated DBA states after the transitions are pruned.

  Simplifying the DBA is important because it can greatly reduce the number of transitions as well as states in the DBA. The impossible transitions can be identified by checking whether their associated propositions will be evaluated false with a given labeling function. For example, if a transition is triggered by the proposition ``a&b`` while the regions labeled ``a`` and ``b`` in the state space of the dynamical system do not intersect (i.e., :math:`L^{-1}(a)\cap L^{-1}(b)=\emptyset`), then this transition is impossible and can be deleted from the DBA. Therefore, trimming transitions are closely related to the definition of the labeling function.

- Extract the metadata of the DBA, i.e., ``nNodes``, ``nAP``, ``acc``, and ``q0``.
- Correctly pass all the DBA transitions to the solvers. As we have seen from how transitions are parsed in ROCS currently, there are two ways to handle this part:

  * Pass all the transitions from the DBA translated by Spot to the variable ``arrayM``, which can be loaded by the ``CSolver`` or ``BSolver``.
  * Write a new function for each solver engine that reads transitions directly from the result translated by Spot (which needs to be simplified first).

.. note::
   This project requires a solid understanding of how the DBA stored and used in Spot. Luckily, Spot is very well documented in this `site <https://spot.lrde.epita.fr/>`_. To integrate the function of Spot, you may find it useful to follow their `C++ code examples <https://spot.lrde.epita.fr/tut.html>`_ and `compiling instructions <https://spot.lrde.epita.fr/compile.html>`_.



Python bindings
---------------
Having a Python interface to the C++ source code in ROCS can facilitate research in formal control synthesis. The user does not need to struggle with writing and compiling a C++ program. We are working on developing Python methods that conveniently

- interpret the DBA translated by Spot Python interface
- define system dynamics
- define labeling functions
- call the two solver engines implemented in C++
- save and display control synthesis results

In this way, writing a C++ `main` program will be replaced by writing a Python script.



Support general LTL formulas
----------------------------
The current limitation of ROCS is that it cannot provide a robustly complete solution for the LTL formula that can only be translated into a non-deterministic Buechi automaton (NBA). Complete control synthesis algorithms rely on a deterministic :math:`\omega` -automaton, and therefore the NBA needs to be determinized to a deterministic Rabin automaton (DRA) or deterministic parity automaton (DPA). To work with DRA or DPA, Rabin or parity game algorithms have to be integrated in order to perform control synthesis. Solving Rabin or parity games are more expensive than solving a Buechi game: the complexity of solving a Rabin game is :math:`\mathcal{O}(n^{2k}k!)`, where :math:`n` is the number of discrete states of the product system and :math:`k` is the number of Rabin pairs; the complexity of solving a parity game is :math:`\mathcal{O}(n^{\log d})`, where :math:`d` is the maximum color.

An abstraction of a nonlinear dynamical system usually contains large numbers of transitions and states. When the given LTL specification is translated into a large DRA or DPA, control synthesis can take an unbearably long time. To solve such control synthesis problems practically, we are working on Rabin and parity control synthesis algorithms that are optimized in both time and space.
