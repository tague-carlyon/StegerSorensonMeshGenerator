from meshStruct import meshStruct
from paramsStruct import paramsStruct

params = paramsStruct.create(foil='NACA0012', gridGenType='Steger-Sorenson', stegerSorenOmega=0.2, debug=False)
mesh = meshStruct.create(params)
mesh.assignInternalConditions()
mesh.linearStretchMeshToWall()
mesh.plotMesh()
mesh.relaxMesh()
mesh.plotMesh()