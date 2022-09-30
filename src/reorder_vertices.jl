

abstract type Singularity end
abstract type Singularity6D <: Singularity end
abstract type Singularity5D <: Singularity end
abstract type Singularity4D <: Singularity end

struct Singularity6DPoint  <: Singularity6D T::SVector{1,Int64}; S::SVector{1,Int64} end
struct Singularity6DEdge   <: Singularity6D T::SVector{2,Int64}; S::SVector{2,Int64} end
struct Singularity6DFace   <: Singularity6D T::SVector{3,Int64}; S::SVector{3,Int64} end
struct Singularity6DVolume <: Singularity6D T::SVector{4,Int64}; S::SVector{4,Int64} end


struct Singularity5DPositiveDistance  <: Singularity5D end
struct Singularity5DPoint  <: Singularity5D T::SVector{1,Int64}; S::SVector{1,Int64} end
struct Singularity5DEdge   <: Singularity5D T::SVector{2,Int64}; S::SVector{2,Int64} end
struct Singularity5DFace   <: Singularity5D T::SVector{3,Int64}; S::SVector{3,Int64} end


struct Singularity4DPoint  <: Singularity4D T::SVector{1,Int64}; S::SVector{1,Int64} end
struct Singularity4DEdge   <: Singularity4D T::SVector{2,Int64}; S::SVector{2,Int64} end
struct Singularity4DFace   <: Singularity4D T::SVector{3,Int64}; S::SVector{3,Int64} end





""""
Tetrahedron-Tetrahedron reordering
"""
#Common Vertex 
function reorder(sing::Singularity6DPoint)

    # Find the permutation P of t and s that make
    # Pt = [P, A1, A2, A3]
    # Ps = [P, B1, B2, B3]
 
    i = sing.T[1]-1
 
    I = circshift([1,2,3,4],-i)
    for k in 1:i 
        reverse!(I,2,4)
    end 
    
    j = sing.S[1]-1
    J = circshift([1,2,3,4],-j)
    for k in 1:j
      
        reverse!(J,2,4)
    end 

    return SVector{4}(I), SVector{4}(J)
end


#Common Edge
function reorder(sing::Singularity6DEdge)

    # Find the permutation P of t and s that make
    # Pt = [P1, A1, A2, P2]
    # Ps = [P1, B1, B2, P2]
  
    I = sort(sing.T)

    I2 = setdiff([1,2,3,4],I)

    for i in 1:I2[2]-I2[1]-1
        reverse!(I2)
    end
    I3 = hcat(I,I2)

    J = sort(sing.S)

    J2 = setdiff([1,2,3,4],J)

    for i in 1:J2[2]-J2[1]-1
        reverse!(J2)
    end
    J3 = hcat(J,J2)

    return SVector{4}(I3), SVector{4}(J3)
end

#Common Face
function reorder(sing::Singularity6DFace)
    # Find the permutation P of t and s that make
    # Pt = [P1, P2, A1, P3]
    # Ps = [P1, P2, A1, P3]

    i = setdiff([1,2,3,4],sing.T)[1]-1
 
    I = circshift([1,2,3,4],-i)
    for k in 1:i
     
        reverse!(I,2,4)
    end 
    I = circshift(I,2)
    
    j = setdiff([1,2,3,4],sing.S)[1]-1
    J = circshift([1,2,3,4],-j)
    for k in 1:j
     
        reverse!(J,2,4)
    end 
    J = circshift(J,2)

    return SVector{4}(I), SVector{4}(J)
end

#Common Volume
function reorder(sing::Singularity6DVolume)

    I = SVector{4,Int64}([1,2,3,4])
    J = SVector{4,Int64}([1,2,3,4])

    return I, J
end

"""
Tetrahedron-Triangle reordering
"""

#Common Volume
function reorder(sing::Singularity5DPositiveDistance)

    I = SVector{4,Int64}([1,2,3,4])
    J = SVector{3,Int64}([1,2,3])

    return I, J
end


#Common Vertex 
function reorder(sing::Singularity5DPoint)

    # Find the permutation P of t and s that make
    # Pt = [P, A1, A2, A3]
    # Ps = [P, B1, B2]
  
    i = sing.T[1]-1
    I = circshift([1,2,3,4],-i)
    for k in 1:i
        reverse!(I,2,4)
    end 
    
    J = circshift([1,2,3],-sing.S[1]+1)

    return SVector{4}(I), SVector{3}(J)
end

#Common Edge
function reorder(sing::Singularity5DEdge)

    # Find the permutation P of t and s that make
    # Pt = [P1, A1, A2, P2]
    # Ps = [P1, B1, P2]
 
    #Tetrahedron
    I = sort(sing.T)

    I2 = setdiff([1,2,3,4],I)

    for i in 1:I2[2]-I2[1]-1
        reverse!(I2)
    end
    I3 = hcat(I,I2)

    #Triangle
    s = setdiff([1,2,3],sing.S)
    J = circshift([1,2,3],-s[1]+2)

    return SVector{4}(I3), SVector{3}(J)
end

#Common Face
function reorder(sing::Singularity5DFace)
    # Find the permutation P of t and s that make
    # Pt = [P1, P2, A1, P3]
    # Ps = [P1, P2, P3]

    #Tetrahedron
    i = setdiff([1,2,3,4],sing.T)[1]-1
 
    I = circshift([1,2,3,4],-i)
    for k in 1:i
     
        reverse!(I,2,4)
    end 
    I = circshift(I,2)
    
    #Triangle
    J = SVector{3,Int64}([1,2,3])

    return SVector{4}(I), J
end


"""
Triangle-Trangle reordering
"""

#Commmon Vertex
function reorder(sing::Singularity4DPoint)

    # Find the permutation P of t and s that make
    # Pt = [P, A1, A2]
    # Ps = [P, B1, B2]
  
    I = circshift([1,2,3],-sing.T[1]+1)
    
    J = circshift([1,2,3],-sing.S[1]+1)

    return I, J
end

#Commmon Edge
function reorder(sing::Singularity4DEdge)

    # Find the permutation P of t and s that make
    # Pt = [P1, A1, P2]
    # Ps = [P1, B1, P2]
  
    t = setdiff([1,2,3],sing.T)
    I = circshift([1,2,3],-t[1]+2)
    s = setdiff([1,2,3],sing.S)
    J = circshift([1,2,3],-s[1]+2)

    return I, J
end

#Commmon Face
function reorder(sing::Singularity4DFace)

    # Find the permutation P of t and s that make
    # Pt = [P1, P2, P3]
    # Ps = [P1, P2, P3]
 
    I = [1,2,3]
    J = [1,2,3]

    return I, J
end



