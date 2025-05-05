from meshStruct import meshStruct
from paramsStruct import paramsStruct

params = paramsStruct.create(foil='NACA0010')
mesh = meshStruct.create(params)
mesh.assignInternalConditions()
mesh.linearStretchMeshToWall()
mesh.plotMesh()
mesh.relaxMesh()
mesh.plotMesh()