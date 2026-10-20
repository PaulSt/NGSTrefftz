from ngsolve import *
from ngsolve.meshes import MakeStructured2DMesh
from xfem import ElementAggregation
from ngstrefftz import AggregatedL2Space, PUFESpace, TrefftzEmbedding, EmbeddedTrefftzFES
from dg import Lap, dgell


def helmholtz_error(fes, test, top, omega):
    exact = exp(1j*(x+y))
    rhs = (2-omega**2)*exact
    embedding = TrefftzEmbedding(top=top, trhs=rhs*test.TestFunction()*dx,
                                 eps=1e-8)
    embedded = EmbeddedTrefftzFES(embedding)
    particular = GridFunction(fes)
    particular.vec.data = embedding.GetParticularSolution()
    a, f = dgell(embedded, exact, rhs, particular, C=-omega**2)
    gfu = GridFunction(embedded)
    gfu.vec.data = a.mat.Inverse(embedded.FreeDofs())*f.vec
    error = gfu+particular-exact
    return sqrt(Integrate(InnerProduct(error, error), fes.mesh)).real


def test_aggregate_helmholtz():
    mesh = MakeStructured2DMesh(nx=3, ny=3, quads=False)
    roots = BitArray([el.nr % 2 == 0 for el in mesh.Elements(VOL)])
    aggregation = ElementAggregation(mesh, roots, ~roots)
    order, omega = 4, 1
    fes = AggregatedL2Space(mesh, aggregation, order=order, complex=True,
                            dgjumps=True)
    test = AggregatedL2Space(mesh, aggregation, order=order-2, complex=True)
    u, v = fes.TrialFunction(), test.TestFunction()
    top = (-Lap(u)-omega**2*u)*v*dx
    assert helmholtz_error(fes, test, top, omega) < 1e-5


def test_pufem_helmholtz():
    mesh = MakeStructured2DMesh(nx=3, ny=3, quads=False)
    order, omega = 4, 1
    fes = PUFESpace(mesh, order=order, complex=True, dgjumps=True)
    test = PUFESpace(mesh, order=order-2, complex=True)
    u, v = fes.TrialFunction(), test.TestFunction()
    top = (-Trace(u.Operator("enrichment_hesse"))
           - omega**2*u.Operator("enrichment"))*v*dx(bonus_intorder=2)
    assert helmholtz_error(fes, test, top, omega) < 1e-5


def test_aggregate_complex_inner_product():
    from embt import _check_inner_product_embedding

    mesh = MakeStructured2DMesh(nx=2, ny=2, quads=False)
    roots = BitArray([el.nr % 2 == 0 for el in mesh.Elements(VOL)])
    aggregation = ElementAggregation(mesh, roots, ~roots)
    fes = AggregatedL2Space(mesh, aggregation, order=4, complex=True)
    test = AggregatedL2Space(mesh, aggregation, order=2, complex=True)
    _check_inner_product_embedding(fes, test, complex_weight=True)


if __name__ == "__main__":
    test_aggregate_helmholtz()
    test_pufem_helmholtz()
    test_aggregate_complex_inner_product()
