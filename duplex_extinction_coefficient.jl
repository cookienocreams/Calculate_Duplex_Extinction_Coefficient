#Calculate the extinction coefficient of an input RNA or DNA duplex

function IsGC(base::Char)
    base == 'G' && return 1
    base == 'C' && return 1
    base == 'A' && return 0
    base == 'T' && return 0
    base == 'U' && return 0
    throw(ErrorException("Unknown nucleotide: $base"))
end

function CalculateDuplexExtinctionCoefficient(First_Strand :: String, Second_Strand :: String
    , First_Strand_ϵ :: Int64, Second_Strand_ϵ :: Int64
    )

    #Need one sequence to find the total percent GC
    combined_sequences = string(First_Strand, Second_Strand)

    duplex_gc :: Float64 = mapreduce(IsGC, +, combined_sequences) / length(combined_sequences)
    duplex_at :: Float64 = 1 - duplex_gc

    hypochromicity :: Float64 = (0.059 * duplex_gc) + (0.287 * duplex_at)

    duplex_extinction_coefficient = (First_Strand_ϵ + Second_Strand_ϵ) * (1 - hypochromicity) |> round

    println(string("The extinction coefficient is ", duplex_extinction_coefficient))
end

function main()
    print("Input the sequence of the first strand in the duplex: ")
    strand_1 :: String = readline()

    print("Input the sequence of the second strand in the duplex: ")
    strand_2 :: String = readline()

    print("Input the extinction coefficient of the first strand in the duplex: ")
    strand_1_ϵ :: Int64 = parse(Int64, readline())

    print("Input the extinction coefficient of the second strand in the duplex: ")
    strand_2_ϵ :: Int64 = parse(Int64, readline())

    CalculateDuplexExtinctionCoefficient(strand_1, strand_2, strand_1_ϵ, strand_2_ϵ)
end

main()