using Thrase
using Test
using Printf


@testset "Thrase.jl" begin
   
   try
      localARGS = ["./examples/test.dat"]
      print(joinpath(@__DIR__,"/src/2D_stripped/stripped_BP1-QD_driver.jl"))
      include("/Users/brittanyerickson/Desktop/Thrase/Thrase.jl/src/2D_stripped/stripped_BP1-QD_driver.jl");
      print("success!\n")
   catch
      print("cannot run bp1-qd variation")
   end

end
