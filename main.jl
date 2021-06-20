# GRAPE ALGORITHM:
# Coder: Caicailiushui
# Date: 20210-05-06
# Description:
# this program is used to generate optimal control signal.

# Input module: pulse shape
# Iteration module: precision control, target function, optimal control

# Output data: density matrix, optimal pulse shape

# include("vars.jl")
include("Tmp.jl")
import .Tmp
Tmp.say_hello()
Hilbert_to_Liouville_index(4)
H = [1 8;
     8 -1]
Tmp.CommutatorMatrix(H)
Tmp.AntiCommutatorMatrix(H)

@show ndim,dt,nslice
VinitArray = zeros(ndim,nslice)

using SpecialFunctions
sigmastd = tmax / (4.0 * sqrt(log(4.0)))
const1 = exp(-tmax^2 / (8.0 *sigmastd^2))
const2 = sqrt(2.0*pi*sigmastd^2) * erf(tmax / (sqrt(8.0)*sigmastd))
const3 = tmax * exp(-tmax^2 / (8.0*sigmastd^2))

ti_tmp = 0.0

for ti = 1:nslice
    ti_tmp = TimeArray[ ti ]
    VinitArray[ 1,ti ] = pi * ( exp(-(ti_tmp + 0.5*dt - 0.5*tmax)^2 / (2.0 * sigmastd^2) ) - const1 ) / (const2 - const3)
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
RhotArray[:,:,1]
PtArray = zeros(ndim,nt)
PtArray[:,1] = [1.0 0 0]

C_operator = zeros(ndim,ndim)
C_operator[ diagind(C_operator) ] = [0 1.0 0]




