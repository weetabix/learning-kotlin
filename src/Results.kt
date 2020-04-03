val rbt = rosalindBiosciToolkit()

val dna = "TGATGTGACGGGTTCGCAGGTGCCCTGTGATTAGAGTAATGAACTTCGTCCTAGCGTAAAAGTGCTCCGGCACATAGAT" +
        "TGCTTTCCCGGGATATATATGCATGTAGGAGATGTTTCCTCATGCTTCCGGTAGGCTAGTCCTACCTCATTAGTCCACTTA" +
        "ACAGTGAACCCTTTCCATCATTAGACGTGGAGAACCATTGTCGCCCCCACAAGCAGGATAAGGAGTTCATCCCCTTAAGCA" +
        "CAAAATCAAATACAAAAGGCAAGTTTGAGACTCTTGTCAAATCTACTCTTTTGCGACTCACAGGCGCAGTTGCAACAATTG" +
        "ATAAGGCTTTGAATCATTAGGAAGACCATAGAGCGAACAAGTCCTCAAACGCGAATGACGGGAGATCACTGCAGACTGTGC" +
        "TAGCTACAGCGCGCCTGCCCGTTCCTCGCCCATGGCCCAAAGAGCCCTAAACCGAGGATGGGTTCGTTTCCGAGGCTCGCG" +
        "TCCTTGCCAGTCGCACCCCTATTAATTGACGATCGTTGGTAAGTTTCCTAAATGCGCATAAGCGACTCTAATAATATGGAC" +
        "GCATCGATAGACACCTAGGTTGTAACCGTGGGGTTAAGACATATGTATTCGACAGGTCGTTTATACTAACATGGTACGTGT" +
        "GCTTCGTGAACTATCTGTCAACTGTAGGTCTCCAACGTGATGAAGTAAGCCCTTAGCGGTGGTTATGTGGCTGGGTTAGTA" +
        "CAGTACGTTTGTGTTTCCGACCCAGTATCTGCTTGCCTATGACGTGCCCCCCAGTCGGAGGAACAAGAATCAAGGCTGCTG" +
        "GATTATCATACGTATGACTGTTAAATCGA"

fun rosalindCountNucleotides() {
    //Problem 1
    println()
    print("Rosalind Count Nucleotides: ")
    println(rbt.countBases(dna))
}

fun rosalindDNACompliment() {
    //Problem 2
    println()
    print("Rosalind Reverse Compliment: ")
    println(rbt.reverseCompliment(dna))
}

fun rosalindGCContent() {
    val strands = rbt.stripFASTA("/home/weetabix/IdeaProjects/OOM Programming Trials/testFileGC-Content")
    println()
    print("Rosalind GC Content: ")
    println(rbt.gcContent(strands))
}

fun rosalindHammingDistance() {
    val pairs = rbt.readHammingPair("/home/weetabix/IdeaProjects/OOM Programming Trials/testFileHamming")
    println()
    print("Rosalind Hamming Distance: ")
    println(rbt.calculateHammingDistance(pairs))
}

fun transcribeDNAToRNA() {
    println()
    print("Rosalind Transcribe DNA to RNA: ")
    println(rbt.transDnaToRna(dna))
}

fun main(args: Array<String>) {
    rosalindCountNucleotides()
    rosalindDNACompliment()
    rosalindGCContent()
    rosalindHammingDistance()
    transcribeDNAToRNA()
}

