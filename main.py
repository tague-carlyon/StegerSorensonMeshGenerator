from meshStruct import meshStruct
from paramsStruct import paramsStruct

params = paramsStruct.create(foil='NACA0010', gridGenType='Steger-Sorenson', stegerSorenOmega=0.02)
mesh = meshStruct.create(params)
mesh.assignInternalConditions()
mesh.linearStretchMeshToWall()
mesh.relaxMesh()
mesh.plotMesh()