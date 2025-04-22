import meshStruct
import paramsStruct


params = paramsStruct.ParamsStruct()
mesh = meshStruct.meshStruct(params)
mesh.assignInternalConditions()
mesh.linearStretchMeshToWall()
mesh.plotMesh()
mesh.relaxMesh()
mesh.plotMesh()