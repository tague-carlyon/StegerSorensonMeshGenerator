from meshStruct import meshStruct
from paramsStruct import paramsStruct

params = paramsStruct.create(foil='NACA6412', gridGenType='Steger-Sorenson', stegerSorenOmega=0.0002, debug=True)
mesh = meshStruct.create(params)
mesh.assignInternalConditions()
mesh.linearStretchMeshToWall()
#mesh.plotMesh()
mesh.relaxMesh()
mesh.plotMesh()