from dolfin import *
from mshr import *
import matplotlib.pyplot as plt

def Fiber(w, h, x, y):
    domain = Rectangle(Point(x - w/2, y - h/2), Point(x + w/2, y + h/2))
    return domain

def Composite(w, h, fw, fh, dx, dy, add=None):
    domain = Rectangle(Point(0, 0), Point(w, h))
    if not add is None:
        domain += add
    nx = int(w / dx)
    ny = int(h / dy)
    for i in range(nx):
        x = w * (i + 0.5) / nx 
        for j in range(ny):
            y = h * (j + 0.5) / ny
            domain.set_subdomain(i * ny + j + 1, Fiber(fw, fh, x, y))
    return domain

def Press(w, h, rx, ry, s):
    domain = Rectangle(Point(-w/2, -h), Point(w/2, 0)) \
           - Ellipse(Point(-w/2, -ry), rx, ry, s) \
           - Ellipse(Point(w/2, -ry), rx, ry, s) \
           - Rectangle(Point(-w/2, -h), Point(-w/2 + rx, -ry)) \
           - Rectangle(Point(w/2 - rx, -h), Point(w/2, -ry))
    return domain
    
def Block(w, h, fw, fh, dx, dy, p, pw):
    presses = CSGTranslation(p, Point(pw/2, 0)) \
            + CSGTranslation(p, Point(w-pw/2, 0)) \
            + CSGTranslation(CSGRotation(p, DOLFIN_PI), Point(w/2, h))
    domain = Composite(w, h, fw, fh, dx, dy, presses)
    return domain

def CoarseBlock(w, h, p, pw):
    presses = CSGTranslation(p, Point(pw/2, 0)) \
            + CSGTranslation(p, Point(w-pw/2, 0)) \
            + CSGTranslation(CSGRotation(p, DOLFIN_PI), Point(w/2, h))
    domain = Rectangle(Point(0, 0), Point(w, h)) + presses
    domain.set_subdomain(1, presses)
    return domain

def SparseBlock(w, h, fw, fh, dx, dy, p, pw, rvew, rveh):
    p = p + Rectangle(Point(-pw/2, 0), Point(pw/2, rveh))
    presses = CSGTranslation(p, Point(pw/2, 0)) \
            + CSGTranslation(p, Point(w-pw/2, 0)) \
            + CSGTranslation(CSGRotation(p, DOLFIN_PI), Point(w/2, h))
    domain = Rectangle(Point(0, 0), Point(w, h)) + presses
    domain.set_subdomain(1, presses)
    nx = int(w / dx)
    ny = int(h / dy)
    for i in range(nx):
        x = w * (i + 0.5) / nx 
        for j in range(ny):
            y = h * (j + 0.5) / ny
            if near(x, pw/2, pw/2) and near(y, 0, rveh) \
            or near(x, w-pw/2, pw/2) and near(y, 0, rveh) \
            or near(x, w/2, pw/2) and near(y, h, rveh):
                domain.set_subdomain(i * ny + j + 2, Fiber(fw, fh, x, y))
    return domain

def save(geo, name, resolution, do_plot=False):
    print(name)
    mesh = generate_mesh(geo, resolution)
    domains = MeshFunction("size_t", mesh, 2, mesh.domains())
    File("meshes/" + name + ".xml") << mesh
    File("meshes/" + name + "_domains.xml") << domains
    File("meshes/" + name + "_domains.pvd") << domains
    if do_plot:
        plt.figure(figsize=(16, 12))
        plot(mesh)
        plot(domains)
    

fiberWidth  = 0.05
fiberHeight = 0.025
densityX    = 0.1
densityY    = 0.1
rveWidth    = 0.1
rveHeight   = 0.1
blockWidth  = 2.0
blockHeight = 0.5
pressWidth  = 0.4
pressHeight = 0.1
roundX      = 0.05
roundY      = 0.05
resolution  = 40
coarse      = 20
segments    = 8

#parameters["lloyd_optimize"] = True
#parameters["mesh_resolution"] = resolution
#parameters["edge_truncate_tolerance"] = 1e-16

rve = Composite(rveWidth, rveHeight, fiberWidth, fiberHeight, densityX, densityY)
save(rve, "rve", sqrt(rveWidth*rveHeight)*resolution, True)

press = Press(pressWidth, pressHeight, roundX, roundY, segments)
save(press, "press", sqrt(pressWidth*pressHeight)*resolution, True)

block = Block(blockWidth, blockHeight, fiberWidth, fiberHeight, densityX, densityX, press, pressWidth)
save(block, "block", sqrt(blockWidth*blockHeight)*resolution, True)

block = CoarseBlock(blockWidth, blockHeight, press, pressWidth)
save(block, "coarse", sqrt(blockWidth*blockHeight)*coarse, True)

sparse = SparseBlock(blockWidth, blockHeight, fiberWidth, fiberHeight, densityX, densityX, press, pressWidth, rveWidth, rveHeight)
save(sparse, "sparse", sqrt(blockWidth*blockHeight)*coarse, True)

