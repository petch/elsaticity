from dolfin import *
import matplotlib.pyplot as plt

def press_bottom(x, on_bound):
    return on_bound and near(x[1], -pressHeight)
 
def press_top(x, on_bound):
    return on_bound and near(x[1], blockHeight + pressHeight)
 
def border_all(x, on_bound):
    return on_bound

class CoefExp(UserExpression):
    def value_shape(s):
        return []
    def eval_cell(s, v, x, c):
        if s.d[c.index] == 0:
            v[0] = s.v1
        else:
            v[0] = s.v2

def Coef(domains, v1, v2):
    W = FunctionSpace(domains.mesh(), "DG", 0)
    C = CoefExp(degree=0)
    C.d = domains
    C.v1 = v1
    C.v2 = v2
    return project(C, W)
    
def epsilon(v):
    return 0.5*(grad(v) + grad(v).T)
 
def sigma(lmd, mu, v):
    return lmd*tr(epsilon(v))*Identity(2)+2*mu*epsilon(v)
            
def problem(name):
    mesh = Mesh("meshes/%s.xml" % name)
    domains = MeshFunction("size_t", mesh, "meshes/%s_domains.xml" % name)
    bounds = MeshFunction("size_t", mesh, 1, 0)
    AutoSubDomain(press_bottom).mark(bounds, 1)
    AutoSubDomain(press_top).mark(bounds, 2)
    
    E   = Coef(domains, E1, E2)
    nu  = Coef(domains, nu1, nu2)
    lmd = nu * E / (1 + nu) / (1 - 2 * nu)
    mu  = E / 2 / (1 + nu)
    
    ds = Measure("ds", subdomain_data=bounds)
    V = VectorFunctionSpace(mesh, "CG", 1)
    u = TrialFunction(V)
    v = TestFunction(V)

    bc = DirichletBC(V, Constant((0, 0)), bounds, 1)
    a = inner(sigma(lmd, mu, u), epsilon(v)) * dx
    L = inner(f, v) * ds(2)
    y = Function(V)
    solve(a == L, y, bc)
    
    File("%s/u.xml" % name) << y
    File("%s/u.pvd" % name) << y
    plt.figure(figsize=(16, 12))
    plot(y)
    return y

def hom_coef(lmd, mu, y, i):
    eps11 =     assemble(epsilon(y)[0, 0]*dx)
    eps12 = 2 * assemble(epsilon(y)[0, 1]*dx)
    eps22 =     assemble(epsilon(y)[1, 1]*dx)
    print("epsilon: ", eps11, eps12, eps22)
    if i == 0:
        epsj = eps11
    elif i == 1:
        epsj = eps12
    else:
        epsj = eps22
    Ej1 = assemble(sigma(lmd, mu, y)[0, 0]*dx)/epsj
    Ej2 = assemble(sigma(lmd, mu, y)[0, 1]*dx)/epsj
    Ej3 = assemble(sigma(lmd, mu, y)[1, 1]*dx)/epsj
    return Ej1, Ej2, Ej3
    
def hom_local(name):
    mesh = Mesh("meshes/%s.xml" % name)
    domains = MeshFunction("size_t", mesh, "meshes/%s_domains.xml" % name)
    
    E   = Coef(domains, E1, E2)
    nu  = Coef(domains, nu1, nu2)
    lmd = nu * E / (1 + nu) / (1 - 2 * nu)
    mu  = E / 2 / (1 + nu)
    
    V = VectorFunctionSpace(mesh, "CG", 1)
    u = TrialFunction(V)
    v = TestFunction(V)

    bcs = [DirichletBC(V, Expression(('x[0]', '0'), degree=1), border_all),
           DirichletBC(V, Expression(('x[1]/2', 'x[0]/2'), degree=1), border_all),
           DirichletBC(V, Expression(('0', 'x[1]'), degree=1), border_all)]

    a = inner(sigma(lmd, mu, u), epsilon(v)) * dx
    L = inner(Constant([0, 0]), v) * ds
    y = Function(V)
    Eh = []
    for i in range(len(bcs)):
        solve(a == L, y, bcs[i])
        File("%s/u%d.xml" % (name, i)) << y
        File("%s/u%d.pvd" % (name, i)) << y
        plt.figure(figsize=(16, 12))
        plot(y)
        Eh.append(hom_coef(lmd, mu, y, i))
    return Eh
   
class ElasticityCoefExp(UserExpression):
    def value_shape(s):
        return (3, 3)
    def eval_cell(s, v, x, c):
        i = s.d[c.index]
        EC = s.EC2
        if i == 0:
            EC = s.EC0
        elif i == 1:
            EC = s.EC1
        for j in range(len(v)):
            v[j] = EC[j]

def ElasticityCoef(domains, EC0, EC1, EC2):
    C = ElasticityCoefExp(degree=0)
    C.d = domains
    C.EC0 = EC0
    C.EC1 = EC1
    C.EC2 = EC2
    return C

def hom_epsilon(v):
    return as_vector((v[0].dx(0), (v[0].dx(1)+v[1].dx(0)), v[1].dx(1)))
    
def hom_coarse(name, Eh):
    mesh = Mesh("meshes/%s.xml" % name)
    domains = MeshFunction("size_t", mesh, "meshes/%s_domains.xml" % name)
    bounds = MeshFunction("size_t", mesh, 1, 0)
    AutoSubDomain(press_bottom).mark(bounds, 1)
    AutoSubDomain(press_top).mark(bounds, 2)
    
    lmd1 = nu1 * E1 / (1 + nu1) / (1 - 2 * nu1)
    mu1  = E1 / 2 / (1 + nu1)
    lmd2 = nu2 * E2 / (1 + nu2) / (1 - 2 * nu2)
    mu2  = E2 / 2 / (1 + nu2)

    EC0 = [Eh[0][0], Eh[1][0], Eh[2][0], Eh[0][1], Eh[1][1], Eh[2][1], Eh[0][2], Eh[1][2], Eh[2][2]]
    EC1 = [lmd1+2*mu1, 0, lmd1, 0, mu1, 0, lmd1, 0, lmd1+2*mu1]
    EC2 = [lmd2+2*mu2, 0, lmd2, 0, mu2, 0, lmd2, 0, lmd2+2*mu2]
    print(EC0)
    print(EC1)
    print(EC2)

    E = ElasticityCoef(domains, EC1, EC1, EC1)
    
    ds = Measure("ds", subdomain_data=bounds)
    V = VectorFunctionSpace(mesh, "CG", 1)
    u = TrialFunction(V)
    v = TestFunction(V)

    bc = DirichletBC(V, Constant((0, 0)), bounds, 1)
    a = inner(E * hom_epsilon(u), hom_epsilon(v)) * dx
    L = inner(f, v) * ds(2)
    y = Function(V)
    solve(a == L, y, bc)
    
    File("%s/u.xml" % name) << y
    File("%s/u.pvd" % name) << y
    plt.figure(figsize=(16, 12))
    plot(y)
    return y
    
blockHeight = 0.5
pressHeight = 0.1
E1, E2 = 4e10, 2e11
nu1, nu2 = 0.15, 0.3
f = Constant((0, -1e5))

#y = problem("block")
Eh = hom_local("rve")
yh = hom_coarse("sparse", Eh)
plt.show()
