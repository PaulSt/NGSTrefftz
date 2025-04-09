.. _documentation: 

Trefftz spaces
========================================
Trefftz finite element spaces extend the ngsolve class of :code:`FESpace`. They provide test and trial functions that are elementwise solutions to chosen PDEs. They are ment to be used in a DG setting with :code:`dgjumps=True`.

.. autoclass:: ngstrefftz.trefftzfespace(mesh,**kwargs)


Embedded Trefftz method
========================================
The embedded Trefftz method produces a Galerkin projection of an underlying discontinuous Galerkin method onto a subspace of Trefftz-type. 
It can be applied to very general cases, including inhomogeneous sources and non-constant coefficient differential operators.
:math:`\newcommand{\Th}{\mathcal{T}_h}
\newcommand{\Vhp}{V^p(\Th)}
\newcommand{\bT}{\mathbf{T}}
\newcommand{\bW}{\mathbf{W}}
\newcommand{\bl}{\mathbf{l}}
\newcommand{\bM}{\mathbf{M}}
\newcommand{\bL}{\mathbf{L}}
\newcommand{\bA}{\mathbf{A}}
\newcommand{\bU}{\mathbf{U}}
\newcommand{\bV}{\mathbf{V}}
\newcommand{\calL}{\mathcal{L}}
\newcommand{\bu}{\mathbf{u}}
\newcommand{\IT}{\mathbb{T}}
\newcommand{\calG}{\mathcal{G}}
\newcommand{\be}{\mathbf{e}}
\newcommand{\bx}{{\mathbf x}}
\newcommand{\inner}[1]{\langle #1 \rangle}
\DeclareMathOperator\Ker{ker}`

Specifically, for an operator :math:`\calL` we are looking to find a basis for its kernel over the space of polynomails. We call this projection :math:`\bT`
Then we can solve the reduced problem: Find :math:`\bu_\IT` so that 

.. math::

    \begin{equation} \label{eq:trefftzlinearsys}
       \bT^T\bA\bT ~ \bu_\IT = \bT^T \bl.
    \end{equation}

The solution in the full polynomial space is then given by :math:`\bu=\bT\bu_\IT`

.. autoclass:: ngstrefftz.TrefftzEmbedding

   .. automethod:: Embed
   .. automethod:: GetEmbedding
   .. automethod:: GetParticularSolution

.. autofunction:: ngstrefftz.EmbeddedTrefftzFES

Trefftz + tent-pitching 
========================================
Given a tent-pitched mesh produced by ngstents the following function returns a solver using Trefftz, or quas-Trefftz, finite elements to solve the acoustic wave equation.

.. autofunction:: ngstrefftz.TWave

The solver class returned by :code:`TWave` is used to set initial and boundary conditions and to propagate the solution in the given tent slab.

.. autoclass:: ngstrefftz.TrefftzTents

   .. automethod:: Propagate
   .. automethod:: SetInitial
   .. automethod:: SetBoundaryCF
