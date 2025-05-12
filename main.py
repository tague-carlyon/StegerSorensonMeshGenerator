from meshStruct import meshStruct
from paramsStruct import paramsStruct

params = paramsStruct.create(foil='NACA6412', gridGenType='Steger-Sorenson', stegerSorenOmega=0.002)
mesh = meshStruct.create(params)
mesh.assignInternalConditions()
mesh.linearStretchMeshToWall()
#mesh.plotMesh()
mesh.relaxMesh()
mesh.plotMesh()