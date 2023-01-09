using Thrase
using Test

my_f(2, 1)

@test my_f(2, 1) == 5

@testset "Thrase.jl" begin

   @test my_f(2, 1) == 5
   @test my_f(6, 1) == 13
  
end
