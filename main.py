import meshStruct
import paramsStruct


params = paramsStruct()
mesh = meshStruct.meshStruct(params)
mesh.assignInternalConditions()
mesh.linearStretchMeshToWall()
mesh.plotMesh()
mesh.relaxMesh()
mesh.plotMesh()