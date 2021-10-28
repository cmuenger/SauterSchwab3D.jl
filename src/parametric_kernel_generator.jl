function generate_integrand_uv(kernel, testref, trialref, testel, trialel)

    function k3(u,v)
        out = @SMatrix zeros(3,3)

        x = neighborhood(testel,u)
        y = neighborhood(trialel,v)

        kernelval = kernel(x,y)
        f = testref(x)
        g = trialref(y)

        jx = jacobian(x)
        jy = jacobian(y)
        ds = jx*jy

        return  SMatrix{3,3}([dot(f[i][1], kernelval*g[j][1])*ds for i=1:3, j=1:3])
    end

    return k3
end