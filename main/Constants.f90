module Constants
    implicit none

    real(8), parameter :: pi = 3.1415926536  
    real(8), parameter :: e = 2.7182818285 
    real(8), parameter :: c = 2.9e10
    complex(16), parameter :: IMG_UNIT = (0, 1) 
    
    real(8), parameter :: R = 3.3
    integer, parameter :: N = 4
    real(8), parameter :: Hr = R / N

    real(8), parameter :: D11 = 10.0
    real(8), parameter :: D22 = 10.0
    real(8), parameter :: D33 = 10.0
    real(8), parameter :: D12 = 10.0

    real(8), parameter :: GParallel = 5.0 * 10e1
    real(8), parameter :: GPerpendicular = 1.0 * 10e2

    real(8), parameter :: Delta1 = 0.0
    real(8), parameter :: Delta2 = 0.0

    real(8), parameter :: HFS = 4.29e10
    real(8), parameter :: q = HFS / c
    real(8), parameter :: V = 1.9825e7 

    real(8), parameter :: Gamma31 = 0.875e7
    real(8), parameter :: Gamma32 = 0.875e7
    real(8), parameter :: GammaSP = Gamma31 + Gamma32
    real(8), parameter :: Gamma = GammaSP + V

    complex(16), parameter :: C1 = (3.0e5, 0.0)
    complex(16), parameter :: C2 = (3.0e5, 0.0)
    real(8), parameter :: a = 1.4

    complex(16), parameter :: DeltaStroke11 = IMG_UNIT * Delta1 + Gamma
    complex(16), parameter :: DeltaStroke12 = -IMG_UNIT * Delta1 + Gamma

    complex(16), parameter :: DeltaStroke21 = IMG_UNIT * Delta2 + Gamma
    complex(16), parameter :: DeltaStroke22 = -IMG_UNIT * Delta2 + Gamma

    complex(16), parameter :: DeltaGamma = IMG_UNIT * (Delta2 - Delta1) + (GParallel + q * q * D12)
    complex(16), parameter :: DeltaGammaStroke = -IMG_UNIT * (Delta2 - Delta1) + (GParallel + q * q * D12)

    ! Dimensions of the matricies used in the project
    integer, parameter :: SYSTEM_VAR_NUM = 8
    integer, parameter :: EXTENDED_MESH_DIM = 8 * (N + 1)
end module Constants