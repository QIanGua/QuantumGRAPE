# GRAPE ALGORITHM:
# Coder: Caicailiushui
# Date: 2021-05-06
# Description:
# this program is used to generate optimal control signal.

# Input module: pulse shape
# Iteration module: precision control, target function, optimal control

# Output data: density matrix, optimal pulse shape

include("vars.jl")
include("Tmp.jl")
import .Tmp
Tmp.say_hello()
Tmp.Hilbert_to_Liouville_index(4)
H = [1 8;
     8 -1];

Tmp.CommutatorMatrix(H)
Tmp.AntiCommutatorMatrix(H)

@show ndim,dt,nslice
VinitArray = zeros(ndim,nslice)

using SpecialFunctions
sigmastd = tmax / (4.0 * sqrt(log(4.0)))
const1 = exp(-tmax^2 / (8.0 *sigmastd^2))
const2 = sqrt(2.0*pi*sigmastd^2) * erf(tmax / (sqrt(8.0)*sigmastd))
const3 = tmax * exp(-tmax^2 / (8.0*sigmastd^2))

for ti = 1:nslice
    ti_tmp = TimeArray[ ti ]
    VinitArray[ 1,ti ] = pi * ( exp(-(ti_tmp + 0.5*dt - 0.5*tmax)^2 / (2.0 * sigmastd^2) ) - const1 ) / (const2 - const3)
    # VinitArray[ 1,ti ] = pi * ( exp(-(TimeArray[ ti ] + 0.5*dt - 0.5*tmax)^2 / (2.0 * sigmastd^2) ) - const1 ) / (const2 - const3)
end

# VinitArray

HtArray = zeros(ndim,ndim,nslice)

# Hd = [0,0; -0.2 * 2.0 * pi]
# diagind(Hd)
# Hd  |> diagind
using LinearAlgebra

Hd = zeros(ndim,ndim)
Hd[diagind(Hd)] = [0.0, 0.0, -0.2 * 2.0 * pi]

SmatArray = zeros(Complex,ndim,ndim,3)

SmatArray[:,:,1] = 0.5*[0 1 0;
                        1 0 sqrt(2.0);
                        0 sqrt(2.0) 0]


SmatArray[:,:,2] = 0.5*[0 -1im 0;
                        1im 0 -1im*sqrt(2.0);
                        0 1im*sqrt(2.0) 0]


SmatArray[:,:,3] = [0 0 0;
                    0 1 0;
                    0 0 2]

Rho_init = [1.0 0 0;
            0 0 0;
            0 0 0]

RhotArray = zeros(ndim,ndim,nt)

RhotArray[:,:,1] = Rho_init
# RhotArray[:,:,1]
PtArray = zeros(ndim,nt)

PtArray[:,1] = [1.0 0 0]

C_operator = zeros(ndim,ndim)
C_operator[ diagind(C_operator) ] = [0 1.0 0]



# iteration to evolute the density matrix
VpreArray = Rho_init
DPhi0Array = zeros(ndim,nslice)
VnextArray = VinitArray

UdtArray = zeros(ndim,ndim,nslice)
UtforwardArray = zeros(ndim,ndim,nt)
UtforwardArray[:,:,1] = I(ndim)

UtbackwardArray = zeros(ndim,ndim,nt)

rhotmax = zeros(ndim,ndim)

# function
loop_step = 0
Phi0 = 0.0
precision
while abs(1 - Phi0) >= precision
    loop_step += 1
    if mod(loop_step,5) == 0 || loop_step < 5 || loop_step == loop_max
        @show loop_step
    end

    for ti in 1:nslice
        HtArray[:,:,ti] = Hd + VnextArray[1,ti]*SmatArray[:,:,1] +
            VnextArray[2,ti]*SmatArray[:,:,2] +
            VnextArray[3,ti]*SmatArray[:,:,3]

        UdtArray[:,:,ti]         = exp(-1im * dt * HtArray[:,:,ti])
        UtforwardArray[:,:,ti+1] = UdtArray[:,:,ti]*UtforwardArray[:,:,ti]
        RhotArray[:,:,ti+1]      = UtforwardArray[:,:,ti+1] * Rho_init * UtforwardArray[:,:,ti+1]

        PtArray[1,ti+1] = RhotArray[1,1,ti+1]
        PtArray[2,ti+1] = RhotArray[2,2,ti+1]
        PtArray[3,ti+1] = RhotArray[3,3,ti+1]
    end

    Ut_max = UtforwardArray[:,:,nt]
    for ti in 1:nt
        UtbackwardArray[:,:,ti] = Ut_max * UtforwardArray[:,:,nt+1-ti]
    end

    rhotmax = RhotArray[:,:,nt]
    Phi0 = real( tr( C_operator' * rhotmax) )
    @show Phi0

    for ti in 1:nslice
        lambda_ti = UtbackwardArray[:,:,nt-ti]' * C_operator * UtbackwardArray[:,:,nt-ti]
        rhoti = RhotArray[:,:,ti + 1]
        for ndim_idx in 1:3
            LHK_commutator = 1im * dt * (SmatArray[:,:,ndim_idx] * rhoti - rhoti * SmatArray[:,:,ndim_idx])
            DPhi0Array[ndim_idx,ti] = - real(tr(lambda_ti' * LHK_commutator))
        end
    end
    VnextArray = VpreArray + epsilon * DPhi0Array
    VpreArray = VnextArray
end


