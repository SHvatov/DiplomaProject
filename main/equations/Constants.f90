module Constants
    implicit none
    real, parameter :: pi = 3.1415926536  
    real, parameter :: e = 2.7182818285 
    complex, parameter :: IMG_UNIT = (0, 1) 

    real, parameter :: q = 0.0
    
    real, parameter :: R = 0.0
    integer, parameter :: N = 1
    real, parameter :: Hr = R / N

    real, parameter :: D11 = 0.0
    real, parameter :: D22 = 0.0
    real, parameter :: D33 = 0.0
    real, parameter :: D12 = 0.0

    real, parameter :: GParallel = 0.0
    real, parameter :: GPerpendicular = 0.0

    real, parameter :: Delta1 = 0.0
    real, parameter :: Delta2 = 0.0

    real, parameter :: HFS = 0.0
    real, parameter :: V = 0.0

    real, parameter :: Gamma31 = 0.0
    real, parameter :: Gamma32 = 0.0
    real, parameter :: GammaSP = Gamma31 + Gamma32
    real, parameter :: Gamma = GammaSP + V

    complex, parameter :: C1 = (0, 0)
    complex, parameter :: C2 = (0, 0)
    real, parameter :: a = 0.0

    complex, parameter :: DeltaStroke11 = IMG_UNIT * Delta1 + Gamma
    complex, parameter :: DeltaStroke12 = -IMG_UNIT * Delta1 + Gamma

    complex, parameter :: DeltaStroke21 = IMG_UNIT * Delta2 + Gamma
    complex, parameter :: DeltaStroke22 = -IMG_UNIT * Delta2 + Gamma

    complex, parameter :: DeltaGamma = IMG_UNIT * (Delta2 - Delta1) + (GParallel + q ** 2 * D12)
    complex, parameter :: DeltaGammaStroke = -IMG_UNIT * (Delta2 - Delta1) + (GParallel + q ** 2 * D12)
end module Constants