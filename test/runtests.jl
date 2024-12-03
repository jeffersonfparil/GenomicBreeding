using GenomicBreeding, Test, Documenter

Documenter.doctest(GenomicBreeding)

@testset "GenomicBreeding.jl" begin
    @test 1 == 0
    genomes::Genomes = Simulate.simulategenomes()
    @test isa(genomes, Genomes)
end
