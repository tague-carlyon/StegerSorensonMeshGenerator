from meshStruct import meshStruct
from paramsStruct import paramsStruct

params = paramsStruct.create(foil='NACA2430', gridGenType='Steger-Sorenson', stegerSorenOmega=0.02, convCriteria=.01,debug=False)
mesh = meshStruct.create(params)
#mesh.plotMesh()
mesh.relaxMesh()
mesh.plotMesh()