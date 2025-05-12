from meshStruct import meshStruct
from paramsStruct import paramsStruct

params = paramsStruct.create(foil='NACA2410', gridGenType='Steger-Sorenson', stegerSorenOmega=0.02)
mesh = meshStruct.create(params)
mesh.assignInternalConditions()
mesh.plotMesh()
mesh.linearStretchMeshToWall()
mesh.plotMesh()
mesh.relaxMesh()
mesh.plotMesh()