module Constants
    implicit none

    real, parameter :: pi = 3.1415926536  
    real, parameter :: e = 2.7182818285 
    real, parameter :: c = 2.9 * 10e10 
    complex, parameter :: IMG_UNIT = (0, 1) 
    
    real, parameter :: R = 3.3
    integer, parameter :: N = 4
    real, parameter :: Hr = R / N

    real, parameter :: D11 = 10.0
    real, parameter :: D22 = 10.0
    real, parameter :: D33 = 10.0
    real, parameter :: D12 = 10.0

    real, parameter :: GParallel = 5.0 * 10e1
    real, parameter :: GPerpendicular = 1.0 * 10e2

    real, parameter :: Delta1 = 0.0
    real, parameter :: Delta2 = 0.0

    real, parameter :: HFS = 4.29e10
    real, parameter :: q = HFS / c
    real, parameter :: V = 1.9825e7 

    real, parameter :: Gamma31 = 0.875e7
    real, parameter :: Gamma32 = 0.875e7
    real, parameter :: GammaSP = Gamma31 + Gamma32
    real, parameter :: Gamma = GammaSP + V

    complex, parameter :: C1 = (3.0e5, 0.0)
    complex, parameter :: C2 = (3.0e5, 0.0)
    real, parameter :: a = 1.4

    complex, parameter :: DeltaStroke11 = IMG_UNIT * Delta1 + Gamma
    complex, parameter :: DeltaStroke12 = -IMG_UNIT * Delta1 + Gamma

    complex, parameter :: DeltaStroke21 = IMG_UNIT * Delta2 + Gamma
    complex, parameter :: DeltaStroke22 = -IMG_UNIT * Delta2 + Gamma

    complex, parameter :: DeltaGamma = IMG_UNIT * (Delta2 - Delta1) + (GParallel + q ** 2 * D12)
    complex, parameter :: DeltaGammaStroke = -IMG_UNIT * (Delta2 - Delta1) + (GParallel + q ** 2 * D12)

    ! Dimensions of the matricies used in the project
    integer, parameter :: SYSTEM_VAR_NUM = 8
    integer, parameter :: EXTENDED_MESH_DIM = 8 * (N + 1)
end module Constants