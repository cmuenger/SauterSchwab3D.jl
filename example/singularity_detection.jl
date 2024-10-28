function singularity_detection(t,s)

    sing = 0;

    D = dimension(t)+dimension(s)
    idx_t = []
    idx_s = []
    for i in 1:length(t)
        v = t[i]
        for j in 1:length(s)
            w = s[j]
            if norm(w-v) < eps(eltype(v)) * 1.0e3
                sing += 1
                push!(idx_t,i)
                push!(idx_s,j)
                break
            end
        end
    end

    if D == 4
        sing == 1 && return SauterSchwab3D.Singularity4DPoint(idx_t,idx_s)
        sing == 2 && return SauterSchwab3D.Singularity4DEdge(idx_t,idx_s)
        sing == 3 && return SauterSchwab3D.Singularity4DFace(idx_t,idx_s)
    elseif D == 5
        sing == 1 && return SauterSchwab3D.Singularity5DPoint(idx_t,idx_s)
        sing == 2 && return SauterSchwab3D.Singularity5DEdge(idx_t,idx_s)
        sing == 3 && return SauterSchwab3D.Singularity5DFace(idx_t,idx_s)
    elseif D == 6
        sing == 1 && return SauterSchwab3D.Singularity6DPoint(idx_t,idx_s)
        sing == 2 && return SauterSchwab3D.Singularity6DEdge(idx_t,idx_s)
        sing == 3 && return SauterSchwab3D.Singularity6DFace(idx_t,idx_s)
        sing == 4 && return SauterSchwab3D.Singularity6DVolume(idx_t,idx_s)
    end
 
end