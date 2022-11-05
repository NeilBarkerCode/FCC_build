# FCC code

def FCCFEA01(X, Y, Z, R, OF, File):
    # X = number of corner spheres in X direction
    # Y = number of corner spheres in Y direction
    # Z = number of corner spheres in Z direction
    # R = Sphere radius
    # OF - overlapping fraction
    # File = location of where spaceclaim file is to be saved
    # Example: File = str('D:/05 Python/02 FCC arrangement/01 spaceclaim python code/')
    # example: FCCFEA01(4, 6, 9, 2, 0.97, File)

    # Packages to be imported
    import numpy as np


    # FCC packing equation
    FCC_Equation = np.sqrt(8) * R * OF
    FCC_Coordination_number = 12
    FCCUnitCellFullSpheres = 4

    # create text file with information in it
    TextFile = open(str(File) + 'FCCFEA01-' + str(X) + '-' + str(Y) + '-' + str(Z) + '-'
                    + str(R) + '-' + str(OF) + '.txt', 'w')

    # create SpaceClaim .py file
    SpaceClaimFile = open(str(File) + 'FCCFEA01-' + str(X) + '-' + str(Y) + '-' + str(Z) + '-'
                          + str(R) + '-' + str(OF) + '.py', 'w')

    # Changes SpaceClaim code to 3D mode
    SpaceClaimFile.write('mode = InteractionMode.Solid' + '\n' + 'result = ViewHelper.SetViewMode(mode, None)' + '\n')

    # create the cube
    SpaceClaimFile.write('result = BlockBody.Create(Point.Create(MM(0), MM(0), MM(0)), Point.Create(MM('
                         + str((X - 1) * FCC_Equation) + '), MM(' + str((Y - 1) * FCC_Equation) + '), MM('
                         + str((Z - 1) * FCC_Equation) + ')), ExtrudeType.ForceAdd, None)' + '\n')

    # Removal of the FCC spheres from the cube
    x1 = np.array([1, 0, 0])
    y1 = np.array([0, 1, 0])
    z1 = np.array([0, 0, 1])

    # corner spheres
    for n1 in range(X):
        for n2 in range(Y):
            for n3 in range(Z):
                xyz1 = n1 * x1 * FCC_Equation + n2 * y1 * FCC_Equation + n3 * z1 * FCC_Equation
                SpaceClaimFile.write('SphereBody.Create(Point.Create(MM(' + str(xyz1[0]) + '), MM(' + str(xyz1[1])
                                     + '), MM(' + str(xyz1[2]) + ')), Point.Create(MM(' + str(xyz1[0] + R)
                                     + '), MM(' + str(xyz1[1]) + '), MM(' + str(xyz1[2])
                                     + ')), ExtrudeType.ForceCut,None)' + '\n')

    # inner sphere X face
    FCC3 = [0 * FCC_Equation, 0.5 * FCC_Equation, 0.5 * FCC_Equation]
    for n10 in range(X):
        for n11 in range(Y - 1):
            for n12 in range(Z - 1):
                FCCNo3 = FCC3 + n10 * x1 * FCC_Equation + n11 * y1 * FCC_Equation + n12 * z1 * FCC_Equation
                SpaceClaimFile.write('SphereBody.Create(Point.Create(MM(' + str(FCCNo3[0]) + '), MM(' + str(FCCNo3[1])
                                     + '), MM(' + str(FCCNo3[2]) + ')), Point.Create(MM(' + str(FCCNo3[0] + R)
                                     + '), MM(' + str(FCCNo3[1]) + '), MM(' + str(FCCNo3[2])
                                     + ')), ExtrudeType.ForceCut,None)' + '\n')

    # inner sphere Y face
    FCC2 = [0.5 * FCC_Equation, 0 * FCC_Equation, 0.5 * FCC_Equation]
    for n7 in range(X - 1):
        for n8 in range(Y):
            for n9 in range(Z - 1):
                FCCNo2 = FCC2 + n7 * x1 * FCC_Equation + n8 * y1 * FCC_Equation + n9 * z1 * FCC_Equation
                SpaceClaimFile.write('SphereBody.Create(Point.Create(MM(' + str(FCCNo2[0]) + '), MM(' + str(FCCNo2[1])
                                     + '), MM(' + str(FCCNo2[2]) + ')), Point.Create(MM(' + str(FCCNo2[0] + R)
                                     + '), MM(' + str(FCCNo2[1]) + '), MM(' + str(FCCNo2[2])
                                     + ')), ExtrudeType.ForceCut,None)' + '\n')

    # inner sphere Z face
    FCC1 = [(0.5 * FCC_Equation), (0.5 * FCC_Equation), (0 * FCC_Equation)]
    # might have to do a separate for loop for each of the internal FCC spheres. add
    for n4 in range(X - 1):
        for n5 in range(Y - 1):
            for n6 in range(Z):
                FCCNo1 = FCC1 + n4 * x1 * FCC_Equation + n5 * y1 * FCC_Equation + n6 * z1 * FCC_Equation
                SpaceClaimFile.write('SphereBody.Create(Point.Create(MM(' + str(FCCNo1[0]) + '), MM('
                                     + str(FCCNo1[1]) + '), MM(' + str(FCCNo1[2]) + ')), Point.Create(MM('
                                     + str(FCCNo1[0] + R) + '), MM(' + str(FCCNo1[1]) + '), MM(' + str(FCCNo1[2])
                                     + ')), ExtrudeType.ForceCut,None)' + '\n')

    # write model information to a text file
    TextFile.write('The Function used: FCCFEA01' + '\n' + '\n')
    TextFile.write('Variables used : ' + '\n' + 'Number of spheres in the X direction ' + str(X)
                   + '\n' + 'Number of spheres in the Y direction ' + str(Y) + '\n'
                   + 'Number of spheres in the Z direction ' + str(Z) + '\n' + 'The sphere radius '
                   + str(R) + 'mm' + '\n' + 'The sphere overlapping fraction ' + str(OF) + '\n' + '\n')

    # calculating the cube dimensions
    XDim = (X - 1) * FCC_Equation
    YDim = (Y - 1) * FCC_Equation
    ZDim = (Z - 1) * FCC_Equation

    TextFile.write('The dimensions of the model are:' + '\n' + 'X = ' + str(XDim) + 'mm' + '\n' + 'Y = ' +
                   str(YDim) + 'mm' + '\n' + 'Z = ' + str(ZDim) + 'mm' + '\n' + '\n')

    # Calculating the unit cell
    VolumeOfSolidUnitCell = FCC_Equation ** 3
    VolumeOfUnitSphere = (4 / 3) * np.pi * R ** 3
    r = R
    D = (2 * R) * OF
    # Window radius of the intersecting spheres
    WindowRadius = float((1 / (2 * D)) * np.sqrt((-D + r - R) * (-D - r + R) * (-D + r + R) * (D + r + R)))
    # finding the volume of the removed lens common to the intersecting spheres
    LensVolume = float((np.pi * (4 * R + D) * ((2 * R - D) ** 2)) / 12)

    # each unit cell has 4 spheres removed and 24 lenz volumes added
    # total number of Lenzs is = coordiantion number / 2 * number of spheres
    VolumeOfUnitCell = VolumeOfSolidUnitCell - (VolumeOfUnitSphere * FCCUnitCellFullSpheres) + (
            LensVolume * FCC_Coordination_number / 2 * FCCUnitCellFullSpheres)

    # Unit volume of the void space

    UnitVoidVolume = VolumeOfSolidUnitCell - VolumeOfUnitCell

    # total volume
    TotalNumberOfUnitCells = (X - 1) * (Y - 1) * (Z - 1)
    TotalSolidVolume = TotalNumberOfUnitCells * VolumeOfSolidUnitCell
    TotalMaterialVolume = TotalNumberOfUnitCells * VolumeOfUnitCell
    TotalVoidVolume = TotalNumberOfUnitCells * UnitVoidVolume
    Porosity = TotalVoidVolume / TotalSolidVolume

    TextFile.write('Total number of unit cells ' + str(TotalNumberOfUnitCells) + '\n')
    TextFile.write('The window radius is ' + str(WindowRadius) + 'mm' + '\n')
    TextFile.write('The solid cube volume is  ' + str(TotalSolidVolume) + 'mm^3' + '\n')
    TextFile.write('The material volume is  ' + str(TotalMaterialVolume) + 'mm^3' + '\n')
    TextFile.write('The void volume is ' + str(TotalVoidVolume) + 'mm^3' + '\n')
    TextFile.write('The volume fraction (porosity) of the cube - ' + str(Porosity) + '\n')

    TextFile.close()

File = str('E:/05 Python/13 testing/')
FCCFEA01(4, 5, 6, 2, 0.97, File)