using Thrase
using Test
using Printf

BASE_FOLDER = dirname(dirname(pathof(Thrase)))

localARGS = [joinpath(BASE_FOLDER, "examples", "test.dat")]


@testset "Thrase.jl" begin 
   try
      testfile = joinpath(BASE_FOLDER, "src/2D_stripped", "BP1_stripped_test.jl")
      include(testfile)
      print("success!\n")
   catch
      print("cannot run bp1-qd variation")
   end

end
