from meshStruct import meshStruct
from paramsStruct import paramsStruct

params = paramsStruct.create(foil='NACA0010', gridGenType='Steger-Sorenson')
mesh = meshStruct.create(params)
mesh.assignInternalConditions()
mesh.linearStretchMeshToWall()
mesh.relaxMesh()
mesh.plotMesh()