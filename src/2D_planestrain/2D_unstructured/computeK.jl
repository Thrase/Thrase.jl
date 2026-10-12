# Two blocks; On an unstructured mesh of quads:

using Thrase
using LinearAlgebra
using OrdinaryDiffEq
using DiffEqCallbacks
using DelimitedFiles
using PyPlot
using SparseArrays



const year_seconds = 31556926
global const ctr = Ref{Int64}(1) 

const RS_FAULT = 7
const VP_FAULT = 8

include("ops_BP3-QD_unstructured.jl")
include("../utils_2D.jl")

const c = .1
const c1 = .1
const b1 = 1
const b2 = 1
const L = 40

#  elastic properties; current set up is homogeneous

function mu(x, z)
        return 2 .+ c .* sin.(b2 .* x ./ L .+ b1 .* z ./ L)
end

function lambda(x, z)
        return 1 .+ c1 .* cos.(b2 .* x ./ L .+ b1 .* z ./ L)
end


function planestrain_solve(Nr, Ns, SBPp, verts, EToV, FToB, EToF, FToE, FToLF, EToO, EToS, nfaces, EToDomain)

    nelems = length(Nr)
    
    ################################## COORDINATE TRANSFORM ###################################
    #
    # Build the local volume operators
    #

    # Create an empty dictionary to store the operators;
    # index via element number, return structure containing local ops.
    OPTYPE = typeof(locoperator(2, 16, 16, mu, lambda))
    lop = Dict{Int64, OPTYPE}()
    
    # Indices corresponding to vertices as you move through blocks/elements
    vstarts = Array{Int64, 1}(undef, nelems + 1)
    vstarts[1] = 1   # start at 1. 

    # Loop over blocks/elements and create local operators:
    for e = 1:nelems

        Np = (Nr[e]+1)*(Ns[e]+1)  # total number of volume points on each element
        vstarts[e+1] = vstarts[e] + Np # fill in vstarts

        # DOMAIN: ordered pairs that define the physical domain
        (x1,x2,x3,x4) = verts[1, EToV[:, e]] 
        (y1,y2,y3,y4) = verts[2, EToV[:, e]]

   
        # Details of transfinite interpolation: straight edges:
        # x-parameterization on each face
        ex = [(α) -> x1 * (1 .- α) / 2 + x3 * (1 .+ α) / 2,
        (α) -> x2 * (1 .- α) / 2 + x4 * (1 .+ α) / 2,
        (α) -> x1 * (1 .- α) / 2 + x2 * (1 .+ α) / 2,
        (α) -> x3 * (1 .- α) / 2 + x4 * (1 .+ α) / 2]

        # derivative of x-parameterization
        exα = [(α) -> -x1 / 2 + x3 / 2,
        (α) -> -x2 / 2 + x4 / 2,
        (α) -> -x1 / 2 + x2 / 2,
        (α) -> -x3 / 2 + x4 / 2]
        
        # y-parameterization on each face
        ey = [(α) -> y1 * (1 .- α) / 2 + y3 * (1 .+ α) / 2,
        (α) -> y2 * (1 .- α) / 2 + y4 * (1 .+ α) / 2,
        (α) -> y1 * (1 .- α) / 2 + y2 * (1 .+ α) / 2,
        (α) -> y3 * (1 .- α) / 2 + y4 * (1 .+ α) / 2]

        # derivative of y-parameterization
        eyα = [(α) -> -y1 / 2 + y3 / 2,
        (α) -> -y2 / 2 + y4 / 2,
        (α) -> -y1 / 2 + y2 / 2,
        (α) -> -y3 / 2 + y4 / 2]

        # Coordinate transformation with transfinite interpolation
        xt(r,s) = transfinite_blend(ex[1], ex[2], ex[3], ex[4], exα[1], exα[2], exα[3], exα[4], r, s)
        zt(r,s) = transfinite_blend(ey[1], ey[2], ey[3], ey[4], eyα[1], eyα[2], eyα[3], eyα[4], r, s)

        # Uncomment if coordinate transform user-specified, for example:
        #xt(r,s) = (2 * r, 2 * ones(size(r)), 0 * ones(size(s)))
        #zt(r,s) = (2 * s, 0 * ones(size(r)), 2 * ones(size(s)))

        metrics = create_metrics(Nr[e], Ns[e], mu, lambda, xt, zt) # create coordinate transform


        ###################################################################### 
        # create local finite difference operators on computational domain:
        lop[e] = locoperator(SBPp, Nr[e], Ns[e], mu, lambda, metrics, FToB[EToF[:, e]]) 
    end

    # Calculate neighboring penalty parameters and add to lop 
    for e = 1:nelems
        lop[e].neighborZ  .= calculate_neighbors(lop, e, FToB, EToF, FToE, FToLF, EToO)
    end
    

    # Assemble the global volume operator and factor:
    Ao = global_operator(lop, vstarts, FToB, FToE, FToLF, EToO, EToS, Nr, Ns)
    A = lu(Ao)

    # construct Tr and BΓ on element 1 (should be the same if computed on element 2):
    Np = (Nr[1]+1)*(Ns[1]+1) # number of volume points on each element 
    N = Ns[1]+1 # number of face points
    
    lf1 = 2;
    lf2 = 1; 

    T1 = lop[1].T[lf1]
    T2 = lop[2].T[lf2]

    # restriction opersators 
    eRST1 = lop[1].eRST[lf1]
    eRST2 = lop[2].eRST[lf2]
    
    sJ1 = lop[1].sJ[lf1]
    sJ2 = lop[2].sJ[lf2]


    Ident = sparse(I, N, N)
    av = [0.5 0.5 0   0;
          0   0   0.5 0.5]
    Avg = kron(av, Ident) # averaging operator
                

    nxf1 = lop[1].nx[lf1]
    nyf1 = lop[1].ny[lf1]
    nxf2 = lop[2].nx[lf2]
    nyf2 = lop[2].ny[lf2]


    spz = spzeros(N, N)
    Normals = [-diagm(nyf1) diagm(nxf1) spz  spz;
          spz  spz -diagm(nyf2) diagm(nxf2);
        -diagm(nxf1) -diagm(nyf1) spz  spz;
          spz  spz -diagm(nxf2) -diagm(nyf2)]


    sPZ = spzeros(N, N^2)
    Trac_op = [ (eRST1 ./ sJ1)*T1[1,1] (eRST1 ./ sJ1)*T1[1,2] sPZ sPZ;
                (eRST1 ./ sJ1)*T1[2,1] (eRST1 ./ sJ1)*T1[2,2] sPZ sPZ;
                sPZ sPZ (eRST2 ./ sJ2)*T2[1,1] (eRST2 ./ sJ2)*T2[1,2];
                sPZ sPZ (eRST2 ./ sJ2)*T2[2,1] (eRST2 ./ sJ2)*T2[2,2]];

    
    Tr = Avg*Normals*Trac_op # Tgamma operator

    
    sJZ11_1 = 0.5*(lop[1].IsJZ[lf1][1, 1])' + (lop[1].neighborZ[lf1][1, 1])'
    sJZ12_1 = 0.5*(lop[1].IsJZ[lf1][1, 2])' + (lop[1].neighborZ[lf1][1, 2])'
    sJZ21_1 = 0.5*(lop[1].IsJZ[lf1][2, 1])' + (lop[1].neighborZ[lf1][2, 1])'
    sJZ22_1 = 0.5*(lop[1].IsJZ[lf1][2, 2])' + (lop[1].neighborZ[lf1][2, 2])'

    sJZ11_2 = 0.5*(lop[2].IsJZ[lf2][1, 1])' + (lop[2].neighborZ[lf2][1, 1])'
    sJZ12_2 = 0.5*(lop[2].IsJZ[lf2][1, 2])' + (lop[2].neighborZ[lf2][1, 2])'
    sJZ21_2 = 0.5*(lop[2].IsJZ[lf2][2, 1])' + (lop[2].neighborZ[lf2][2, 1])'
    sJZ22_2 = 0.5*(lop[2].IsJZ[lf2][2, 2])' + (lop[2].neighborZ[lf2][2, 2])'

  
    A11_1 =  -lop[1].Hinv * sJZ11_1 * lop[1].eRS[lf1] * lop[1].H[lf1] + 0.5 * lop[1].Hinv * lop[1].TT[lf1][1, 1] * lop[1].eRS[lf1] * lop[1].H[lf1]
    A21_1 =  -lop[1].Hinv * sJZ21_1 * lop[1].eRS[lf1] * lop[1].H[lf1] + 0.5 * lop[1].Hinv * lop[1].TT[lf1][1, 2] * lop[1].eRS[lf1] * lop[1].H[lf1]
    A12_1 =  -lop[1].Hinv * sJZ12_1 * lop[1].eRS[lf1] * lop[1].H[lf1] + 0.5 * lop[1].Hinv * lop[1].TT[lf1][2, 1] * lop[1].eRS[lf1] * lop[1].H[lf1]
    A22_1 =  -lop[1].Hinv * sJZ22_1 * lop[1].eRS[lf1] * lop[1].H[lf1] + 0.5 * lop[1].Hinv * lop[1].TT[lf1][2, 2] * lop[1].eRS[lf1] * lop[1].H[lf1]
      
    A11_2 =  -lop[2].Hinv * sJZ11_2 * lop[2].eRS[lf2] * lop[2].H[lf2] + 0.5 * lop[2].Hinv * lop[2].TT[lf2][1, 1] * lop[2].eRS[lf2] * lop[2].H[lf2]
    A21_2 =  -lop[2].Hinv * sJZ21_2 * lop[2].eRS[lf2] * lop[2].H[lf2] + 0.5 * lop[2].Hinv * lop[2].TT[lf2][1, 2] * lop[2].eRS[lf2] * lop[2].H[lf2]
    A12_2 =  -lop[2].Hinv * sJZ12_2 * lop[2].eRS[lf2] * lop[2].H[lf2] + 0.5 * lop[2].Hinv * lop[2].TT[lf2][2, 1] * lop[2].eRS[lf2] * lop[2].H[lf2]
    A22_2 =  -lop[2].Hinv * sJZ22_2 * lop[2].eRS[lf2] * lop[2].H[lf2] + 0.5 * lop[2].Hinv * lop[2].TT[lf2][2, 2] * lop[2].eRS[lf2] * lop[2].H[lf2]

    BΓ =   [A11_1 A21_1;
            A12_1 A22_1;
            A11_2 A21_2;
            A12_2 A22_2]
    
    K = Tr * (A \ Matrix(BΓ))


    return K

end


let 
    # SBP interior order of accuracy
    SBPp = 2


    # Set up geometry (2 blocks)
    verts = [-40 0 40 -40 0 40;
              -80 -80 -80 0 0 0];

    FToB = [1; 7; 1; 2; 2; 2; 2];

     EToV = [1 2; 
             2 3; 
             4 5;
             5 6];
     EToF = [1 2; 
             2 3;
             4 5;
             6 7];
    
   EToDomain = [1 2]


    # number of elements and faces
    (nelems, nfaces) = (size(EToV, 2), size(FToB, 1))
    (FToE, FToLF, EToO, EToS) = connectivityarrays(EToV, EToF)

    NN = [11 41 201] # number of nodes in each direction in each block. 

    for n in NN

        Nr = [n n]
        Ns = [n n]

        K = planestrain_solve(Nr, Ns, SBPp, verts, EToV, FToB, EToF, FToE, FToLF, EToO, EToS, nfaces, EToDomain) 

        σ = svdvals(K)

        plt.scatter(1:length(σ), log.(σ))
        
        @show σ
        
    end
end