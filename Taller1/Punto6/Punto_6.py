# Scene
from yade import plot
from yade import qt
yade.qt.Renderer()
qt.View()

# Material
O.materials.append(FrictMat(young=62e9, poisson=0.2, density=2230, label='reloj'))
O.materials.append(FrictMat(young=75e9, poisson=0.16, density=2634, label='arena'))


# Import a mesh
from yade import ymport

# Can also add .stl, .geo and more files
id_HourGl = O.bodies.append(ymport.gmsh("r_0=1.5.mesh", scale=1000.0, color=(0, 0, 1),material='reloj'))



from yade import pack

sp = pack.SpherePack()


sp.makeCloud((-47, -47, 95), (59, 59, 194), num = 400, rMean=5.5, rRelFuzz=0, seed=2)

sp.toSimulation(material='arena')

#Create a Wall
Wall = yade.utils.wall(Vector3(0,0,0), 2, sense=1, material='reloj')
Wall_ID  = O.bodies.append(Wall)


# Create the engine

O.dt = 1e-3

O.engines = [
    ForceResetter(),
    InsertionSortCollider([Bo1_Sphere_Aabb(), Bo1_Facet_Aabb(), Bo1_Wall_Aabb()]), # Contact Detection
    InteractionLoop(
            [Ig2_Sphere_Sphere_ScGeom(), Ig2_Facet_Sphere_ScGeom(), Ig2_Wall_Sphere_ScGeom()],
            [Ip2_FrictMat_FrictMat_MindlinPhys(en=0.07)],
            [Law2_ScGeom_MindlinPhys_Mindlin()],
            
    ),
    NewtonIntegrator(damping=0, gravity=[0, 0, -10.0]),
]

def piso():
    if O.iter > 6000:
        # Se cambia la posicion de la pared
        O.bodies[Wall_ID].state.pos=Vector3(0,0,200)

O.engines += [PyRunner(command='piso()', iterPeriod=1)]


prev_cont = 0 
def flujo():
    global prev_cont
    if O.iter > 6000:
        cont = 0
         
        for body in O.bodies:
            z = body.state.pos[2]
            if isinstance(body.shape, yade.utils.Sphere) and z < 0:
                cont += 1
        if cont > prev_cont:
            with open('r_1.5--d_5.5.txt', 'a', encoding='utf-8') as file:
                file.write(f"{O.iter}, {cont}\n")
            prev_cont = cont


O.engines += [PyRunner(command='flujo()', iterPeriod=1)]


O.saveTmp()
