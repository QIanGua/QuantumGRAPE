using Parameters

@consts begin
    tmin = 0.0
    tmax = 10.0
    dt   = 0.5
    ndim = 3
    epsilon = 0.1
    precision = 0.01
    loop_max = 50
end


@with_kw struct vars
    TimeArray = collect(tmin:dt:tmax)
    nt = length(TimeArray)
    nslice = nt - 1
end

TimeArray = vars().TimeArray
nt = vars().nt
nslice = vars().nslice

