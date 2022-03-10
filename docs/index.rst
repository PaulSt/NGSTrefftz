.. NGSTrefftz documentation master file, created by
   sphinx-quickstart on Wed Mar  2 15:24:26 2022.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. mdinclude:: ../README.md

.. .. include:: notebooks/laplace.rst


Introduction to Trefftz-DG
========================================


.. jupyter-execute::
    :hide-code:
    :hide-output:

    from ngsolve import *
    from ngstrefftz import *
    from netgen.geom2d import unit_square

.. jupyter-execute::
    :hide-output:

    mesh = Mesh(unit_square.GenerateMesh(maxh=.3))
    fes = trefftzfespace(mesh,order=4,eq="laplace")

.. jupyter-execute::
    :hide-code:
    :hide-output:

    exact = exp(x)*sin(y)
    bndc = exact

    def lap_problem(fes):
        mesh = fes.mesh
        order = fes.globalorder
        alpha = 4
        n = specialcf.normal(mesh.dim)
        h = specialcf.mesh_size
        u = fes.TrialFunction()
        v = fes.TestFunction()

        jump_u = (u-u.Other())*n
        jump_v = (v-v.Other())*n
        mean_dudn = 0.5 * (grad(u)+grad(u.Other()))
        mean_dvdn = 0.5 * (grad(v)+grad(v.Other()))

        a = BilinearForm(fes,symmetric=True)
        a += grad(u)*grad(v) * dx \
            +alpha*order**2/h*jump_u*jump_v * dx(skeleton=True) \
            +(-mean_dudn*jump_v-mean_dvdn*jump_u) * dx(skeleton=True) \
            +alpha*order**2/h*u*v * ds(skeleton=True) \
            +(-n*grad(u)*v-n*grad(v)*u)* ds(skeleton=True)

        f = LinearForm(fes)
        f += alpha*order**2/h*bndc*v * ds(skeleton=True) \
             +(-n*grad(v)*bndc)* ds(skeleton=True)

        with TaskManager():
            a.Assemble()
            f.Assemble()
        return (a,f)

.. jupyter-execute::

    errors = []
    for p in range(1,6):
        errors.append([])
        mesh = Mesh(unit_square.GenerateMesh(maxh=1))
        for h in range(6):
            fes = trefftzfespace(mesh,order=p,eq="laplace")
            (a,f) = lap_problem(fes)
            gfu = GridFunction(fes)
            with TaskManager():
                gfu.vec.data = a.mat.Inverse(inverse='sparsecholesky') * f.vec
            errors[p-1].append(sqrt(Integrate((gfu-exact)**2, mesh)))
            mesh.Refine()
    print(errors)


.. jupyter-execute::
    :hide-code:

    import matplotlib.pyplot as plt

    # plt.subplots_adjust(hspace=1)

    # log y axis
    plt.subplot(121)
    plt.semilogy()
    plt.title('semilogy')
    plt.grid(True)

    # log x and y axis
    plt.subplot(122)
    for p in range(1,6):
        plt.loglog([0.5**i for i in range(6)], errors[p-1],'-o',label=f'p = {p}')
    ax=plt.gca()
    ax.invert_xaxis()
    plt.grid(True)
    plt.legend()
    plt.title('loglog base 2 on x')

    plt.show()



Documentation
========================================

Trefftz spaces
----------------------------------------
Trefftz finite element spaces extend the ngsolve class of :code:`FESpace`. They provide test and trial functions that are elementwise solutions to chosen PDEs. They are ment to be used in a DG setting with :code:`dgjumps=True`.

.. autoclass:: ngstrefftz.trefftzfespace(mesh,**kwargs)


Embedded Trefftz method
----------------------------------------
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

Specifically, for an operator :math:`\calL` we are looking to find a basis for its kernel over the space of polynomails. We call this projectoin :math:`\bT`
Then we can solve the reduced problem: Find :math:`\bu_\IT` so that 

.. math::

    \begin{equation} \label{eq:trefftzlinearsys}
       \bT^T\bA\bT ~ \bu_\IT = \bT^T \bl.
    \end{equation}

The solution in the full polynomial space is then given by :math:`\bu=\bT\bu_\IT`


.. autofunction:: ngstrefftz.TrefftzEmbedding

Trefftz + tent-pitching 
----------------------------------------
Given a tent-pitched mesh produced by ngstrefftz the following function returns a solver using Trefftz, or quas-Trefftz, finite elements to solve the acoustic wave equation.

.. autofunction:: ngstrefftz.TWave

The solver class returned by :code:`TWave` is used to set initial and boundary conditions and to propagate the solution in the given tent slab.

.. autoclass:: ngstrefftz.TrefftzTents

   .. automethod:: Propagate
   .. automethod:: SetInitial
   .. automethod:: SetBoundaryCF


NGSTrefftz Notebooks
========================================

.. toctree::
   :maxdepth: 1

   notebooks/index.ipynb
   notebooks/laplace.ipynb
   notebooks/helmholtz.ipynb
   notebooks/twave.ipynb
   notebooks/qtwave.ipynb
   notebooks/twavetents.ipynb
   notebooks/embTrefftz.ipynb
   notebooks/embTrefftz-poi.ipynb
   notebooks/embTrefftz-wave.ipynb
   notebooks/embTrefftz-helm.ipynb
   notebooks/embTrefftz-adv.ipynb
