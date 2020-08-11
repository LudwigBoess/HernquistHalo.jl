module HernquistHalo

    using ProgressMeter
    using Roots

    export HernquistParameters, 
            sample_halo

    struct HernquistParameters
        Npart::Int64
        Mtot::Float64
        r0::Float64
        rmax::Float64
    end

    global const Hq_const = 1.20223157581242

    function sample_xr_xv(Mtot::Float64, r0::Float64, rmax::Float64)

        xr = 0.0
        xv = 0.0
        e = 10.0
        while ( e > 0.0 || e < (-Mtot/r0) )

            # sample r from system size
            xr = rmax * rand()
            # sample v from escape velocity
            xv = sqrt(2Mtot/r0) * rand()

            # Make sure energy is negative -> particle is bound!
            e = 0.5xv^2 - Mtot/(xr+r0)
        end

        return xr, xv, e
    end

    function get_pn(xr::Float64, xv::Float64, e::Float64, Mtot::Float64, r0::Float64)

        q = sqrt(-1.0*r0*e/Mtot)
        q2 = q*q

        fq = ( 3.0*asin(q) + q*sqrt(1.0 - q2) *
             ( 1.0 - 2.0 * q2) * (8.0*q2*q2 - 8.0*q2 - 3.0)) /
             ( 1.0 - q2 )^2.5

        return ( xr*xr/(r0*r0) ) * ( xv*xv / ( Mtot/r0 )) * fq / Hq_const
    end

    function sample_particle(xr::Float64, xv::Float64)

        # sample position
        cth   = 2.0 * ( rand() - 0.5 )
        signs = 2.0 * ( rand() - 0.5 )
        sth   = signs/abs(signs) * sqrt( 1.0 - cth*cth)
        phi   = 2π * rand()

        x = xr*sth*cos(phi)
        y = xr*sth*sin(phi)
        z = xr*cth

        # sample velocity
        cth   = 2.0 * ( rand() - 0.5 )
        signs = 2.0 * ( rand() - 0.5 )
        sth   = signs/abs(signs) * sqrt( 1.0 - cth*cth)
        phi   = 2π * rand()

        vx = xv*sth*cos(phi)
        vy = xv*sth*sin(phi)
        vz = xv*cth

        return x, y, z, vx, vy, vz
    end

    function sample_hernquist_pos_vel(Npart::Int64, Mtot::Float64, r0::Float64, rmax::Float64)

        # calculate particle mass
        mk = Mtot*rmax^2 / (rmax + r0) / Npart
        m = mk .* ones(Npart)

        x  = zeros(Npart)
        y  = zeros(Npart)
        z  = zeros(Npart)

        vx = zeros(Npart)
        vy = zeros(Npart)
        vz = zeros(Npart)

        P = Progress(Npart)
        idx = 0

        Np = 1
        while Np < Npart

            xr, xv, e = sample_xr_xv(Mtot, r0, rmax)

            if (e > 0.0)
                error("e = $e")
            end

            pn = get_pn(xr, xv, e, Mtot, r0)

            if (rand() <= pn)
                x[Np], y[Np], z[Np], vx[Np], vy[Np], vz[Np] = sample_particle(xr, xv)

                Np  += 1
                idx += 1
                ProgressMeter.update!(P, idx)

            end
        end

        return [x y z], [vx vy vz]
    end

    function calculate_hernquist_masses(Npart::Int64, Mtot::Float64, r0::Float64, rmax::Float64)
        mk = Mtot * rmax * rmax / (rmax + r0)^2 / Npart
        return mk .* ones(Npart)
    end

    function sample_halo(halo::HernquistParameters)

        m = calculate_hernquist_masses(halo.Npart, halo.Mtot, halo.r0, halo.rmax)

        pos, vel = sample_hernquist_pos_vel(halo.Npart, halo.Mtot, halo.r0, halo.rmax)

        return pos, vel, m
    end



end # module
