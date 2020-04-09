import com.marcinmoskala.math.permutations
import com.marcinmoskala.math.permutationsNumber
//import java.net.URL

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

val text =
    "TGATACTATTATACTATGATACTATTAGATACTATCCGGTATATACTATACATACTATTCATACTATCGATACTATGCGGACCACCGATACTATTGGCTATACTATATACTATATACTATATACTATCGAATACTATGCTCATACTATCTTACAATGGATGAAGCTCAGGACATACTATATGCATACTATGGATACTATGATACTATAGACGATACTATATGGATACTATTTTCATACTATTTTGTCATACTATCTATACTATGTCGTTATACTATGTGAGATACTATCATACTATGGAGGATACTATCATACTATATACTATGCAATACTATATACTATGAATACTATATACTATAATACTATCATACTATGATACTATTAAGGAACATACTATTATACTATTGCGCATACTATATACTATATACTATGGGTCGATACTATATACTATCATACTATATACTATAAGATACTATATACTATTTCAGGATACTATTATACTATATACTATCACCATACTATAATACTATGCGATACTATCAATACTATAGAATACTATATACTATATACTATAATACTATGATACTATATACTATCCGGGCGAATACTATATACTATGGGGTGGATACTATATACTATGTGTAATACTATCAATACTATATTCCATACTATTGGGCTTTCAAATACTATAAATACTATGGATACTATGGAATACTATAATACTATTATACTATCGTGATACTATCATCGTATACTATATACTATGAATACTATGATACTATAATACTATGGTATACTATTGCCAAATACTATATACTATATACTATGCCACGCGCTGAATACTATAATACTATTAGATACTATATACTAT"
val subs = "ATACTATAT"
val rna2protein =
    "AUGACCCCUAACAGUCACACUGCGGCUCGAAUCUCAAUCGUGUGGGGAAUGCUCACUCAUGCGAUCUGUUUUCGAUUGCCCAGACUACAGUGCAAUGACUUCCCCACUGGUUGGGACGAGCUCGCUAAGAAAGGAACAGCAGCAGCGCAGCUAAUACUACUCGCCUUAGAAAAAAGCAGGGAGUACGAGCUCUUUCUACCAUGCGUGCUGCUUCAGUGGCCAUCCGAGUACGUGCGCGUAAAGCGGCGUGGGAGUCCAAGUGCGCAUAGGGGACCCACCGUUAGCAAUAUUGUCUGGUGGGGUAUCAUUAUAUGCGCCUCUAAGGAGGGCCCCUCUGGCGCUCCAACUGAUCCAGGUCGUCGUCAUAUAUAUCGUGUCGAAGCUAAUCAAGCCCUGUUAACUUGCAUCUGUCUUCUAAAGAUACAAACUGUGUCAACUCGCACUGAUUGUCACUUGGCAAGCGAUGACUACAAGCUCCGCAUCGCCGCACUGACGAAAUCCAUGAAACUCACACUGGAUAUGUUCGUCUGGUGGAUCCCCCUCCCCUUGAGGUCAAUAGACCGAACUACUUUCAGAUUAAGCAGGGGGCAAGGCUGCUCCCCCGUUGCUACGAAAGCCUACGUACGGCCCAGUGAGGGCCCUCGGUGUCCACAGAGUGAUGGUGGUAUAUGUUCGUCAUCGCGAAUUGCGGUUGGCGUCUGCCUCAGCCAGUGUUGGUCAAACAAGCGCGCUUCCAGACGUGACAUUCAGCGGAACUACACAUCUCUGUUGACAAAGACCACCUUGGUAUUGCCGCCGUCCGGAGUCGAGGUUGAGCGGCCCUGGUGCAGUUCAUGCAGAUCCUACGUAGUUCCCACAGUACUCAAGGUUCACCUACAUUUAACGCUAGUUUCAGCUACGAGAGGACCCGACCUCCGCCUCGGGUUACCGCUUGGUGGUGUCCACUGUAGGUUAUAUUAUCGACAAUGUAAUGAAAUUUCUCGUUUACAGAGCCUUUGCCGCAGGGCGACAGUAUUUCGCACGAGGAUCAUAAGCCAUUCGAUGUGUAACUUCGUCUACCGCUCACGAAAAAUGACGUGGAGAGCUAACCAGUCACGGGAUAACCCUACUCAGCGGCGAAGGCCCUGCCCACGGCUUAAAGAAGGGACUCUACUCCGCUGUGUGACUAUGUCAAGGGCUUGGCCAGCCUUGUUGCCAGGCAUAUUAACAGGUCUCUCUAUUAACGCCUACUGCCCUGUUCGGGACUGCCUUUUUGCUACAGCUGGUUUUUACUCUAGAGAAUCCAUUUUCGAGAUAGUCGCGGGCGACCCGCAUACCAUACACGAGGAAGUUCUCCAAAGCGGGUGUACCAUAUUUCCUCGGGCUCCAUUUGAUUACAAACGAGGCCUAACUUCUUCGUUACCUCCAGUUUCAAUGGAGGCCUGGCGGCUUCACAACACACGUUGGGUCACCCCGGAGUGUGCUAUUGGACCAGAGGUCUACGCCACUUGGGGGCCUUCUCAUGCGAUCUCUCCGCAGGUACCCGGGCUGCAGCGACGUGCAACCCCUUUCGUUACACACUUGGACCGUUUAAACGCUGCAGGCCACAGUAUUACAAUGCUAAUUCGGUUGGUUGGGAGAGCGUGCGUGGACAUGUCCCGGCUGGAUUCAGCAUCUCCGGACGCUUUGGGCCAAUGUACCAUUCUCAUGUGCGCCAGGGUUUCCAGGAGUUAUAACUACCUUUUGUGUUAUUGGAAGGCCUUUUCUUGCGUGGAAUCGCAUUUGCGGAGCACGGGGUUGGCCGAGCUGAUCCGCCCUGGUACGGCGACAAAAAUACGGUGGCUGAAGCUGGACCACCAUGAACAGAUAUCAGUAGAGGUCCGCGUUAGGCCACAUCCCUUCGUGAGAUUAUCAAUGGCCAUACUGGGCAUCGAGGCCCCAACUCGCAGAUACCGUGAAAAGGUACCGACUACACAUGCUGGACAUCGUAUACGCAGUAGAAGAGGUUUCGAACACACGGGGCGUACGGGAGAGAAUAGAAGCGCUUGGGUUUGUGCGAAAAGCUGUUUAUGGCUCCAGGGCGGAACAGGAGCAGGGGAAACAGCACUAAUAGCCGCGCCUUUGGAGCACAGCUCGUCACCAAGUAGAGAAAACGUCGAUUGCUUUUCCAAUUGUACAGUCAUGAAGACGCGGCUAUAUAGUGGAUCGAGUACUUUUUGGGUGUACCCGAUUACAACGAUAAUACUGAAAGACACAGGGACCCGUAUUAUCCCCUCCAGCCUAAGCCUUUCAGACAGGAGAGGUAAGCUAGAACAAGGUCAGAAAGAUUACAUGUUGGAUUACAGCAUGUUGAUACGUUUUUCAGAUCCUCAAGAACGUGCCUCCGGGCCCUCGAUGGUAGGGUGUCUCUCGGGCACCUCCAAUACAGCCCAUCGAGUUGGUUGGACGUUAUUUGGUGUGAUUAAGUUGCUUUCCAAGAAACAACGAUUCAAGUCUAUCAAAUUGUCCUCCCGAACAUGGCACGGCCUCCAUACUGGCAUUGGAAGCUGGUACCGCACCGCAAUGAGUACCCCUGGAUGUAUUAUAUCCUUGAACCCCAUCCUGACGCACUCCGUUAAUCAUACGUCGCGAUUAGAGGUGGUCUUACUCGACGCCAUGCGCAUCACACUGGACUCAGAGCAGCAAACGAUCCUGGAUACCCAAUUGAGGAUAAAAGGUCGACUAGCUACUGGGUUGCUGAUAUUGAUUUUUAAGUUAUCUCGGGUAUCAAGACCCGCCCGUGUAACAAUAAGAUGCACACUAGGCCAUCGCCUCGGGGUUCGUGCUAAGUUGUAUGGCGACAUUGUCCAGCGAUGUUGGCACUCCCGCCUGCCUUAUAUAGCCACCUUAGACAAAGGGCAUCACUCGCAAGUUACUUACAUGAAACGCCACUCCGCCUGUCGGCGAAGGCUGCAUUCCACUACCAGGACGGCUACGGUAUACAUGCAGGAUCUGUGGCACCAACGCGGGUCAAGUUGGACGUACAAAAGGUUGGCCGGCUUGAUAGAGGAUCCAAGAUUUGCGGAUUUGGCCUACGAACGCGGGCGCCGUAUGCGGCGGUAUAUCAUGAACGGUCUUCAUCUGAUACCCCAAGCGACAAUAGUUUCUGCACGGCAAGCUGCGUUAGGAAUUGAACCACCCACGUGUCCCCGGGUACGACGCUCCGGAGUGCUAGUUGAACGUCCCAGGUCUGUGCUAAGGAUAGCUCAACAGUGGAUACACAAUGUAGGGCGAUACGUUCGAGAGUCCACUGGAGCAGCUGACCAUAUGCUUAUCCAGAGACGGGAACAUGUUGAUCACUUGGAAGCAACAGACGAGAGAUUAUCGCGCCACAAAAAGAAGACAUGUAAGUUACGAGAGAGAGUGGCCGUAGUAGGUCACGUGCACAGCUCAAUCGGAGGAACUGGACGGAGCUCCAUGCCACAAGUUCCGGUUCAAGGUCCCUUUGGCACUAGCCUCAAUUGUCAGACCACGAACUGUCGCGUCUGGUCUUGGUCCGUUUCGGACUUAUAUUCAGUGGAGCCUCGUGCCUCACAGACCAGUUUCCCCUACGCCAGGAAGGCCCAGACAACCCAGCGCCGAGCCUGCCUGAAUAGGGGUCUAUUUGGUCGGCCUACCCUCUACGGCACUAUCUUAAAGAAGGGAUAUGUAGUUAACGGCGCGCGGCGCUAUGAGAACAUGAGGAUCAGCUGCUCGGACGAUUCUGGGAUGGGGGUACGAGUCAAGGCACGAUCUUGUGCAAGUGCGAUUAUCGUUUUACAACGGUGGUGGUUAGGUGGCGGAUUAACCGGCGAUCAUGCAAAACGAAUCGUGUCGGCCCCCAUGCCGCGAUGGUCGGAUUUUCACGUCGCUUUCGACCCACUAUUCCCGUGUUUUCGUCAUGACCUGAACCCACAUCAUAUCUCUCUUCUGGAGGGUCACGGUAUUCCGGUACAUCAAACCGGAACAUACCCAGAAAUCUUCGGGUACGCAGGCAACAGUAGCAGAGUUAUUACGCACGAUAAAUCAGCCAUUGCCGACAACCUUCACACCAUGUUUAAAUUGGGAGAUAGACGCGGGGCAGGGUCCAAAUCUAUACUCAAAUGUUUAUGGCGCAUCACAAGCCUGAGUCCGAAUGUCUGUUGCGUUCUGAGACGAAGAAAGUGCAAACAAAGCUCUGGGGUCACCCGUAGCCCCCGGUUAUUGGCGUCAGAUAGCACUGUAUCUAGGAUCAUACGGUGGGAAUGUGUAGUACCUCGAGUUUCUAUUCGUGGUCCCCCCCCACAAACGUGGAUCUUGCGUAUAACGACAGGUAUUGGGGUAGUGUGUUUGUGGCGUAUCUCGCAACUACUACCGACUGAGGACCUAAUGUCUUUCUUUAACCGCGGCCGAAACUCUGUCACAGCCAUCCUCAUUUUCCACAGAAGCUGUGUCACGCACAGCACCUUCGAGCAUGAAAGCAAGUACUUUCGCCCGCUUCACCAUGGAAGUCGUGGAUCAUCGAUUCGGAGAAUCGGACAUCCAGGAUCCAACAACGGGGUCGAUAUUUUUGCAUCGGAUAUCCAACCGAUAUGGCAGUGGCGUCAAAUAGCAGCGCCAGCCACUCUCGAGGAGUCAGGACAUAAGAUCAAUCUCGUUUUUAGCAUGCCCAGAUCGUUCAAUCAAACAUCGGAACCGUGGCGGCUUCAUAUGCCGUCUGCCCUGCGCAUGCCCUACGGUCACAGGAUAACUCGGAAUCUGGAAAUUAUCACUGCGUUGGACAUACGGGAUAAUUUCUCGCUUCAAUUCGGCAGCGGGACGGUUCGAAUACACCAUAUUUGCCCAGCUCGUGGGUUUCGUUCGAUAGGUAGGAAGUGGGCUUGGUGGGGGCUUAUGUCAGCGAUUUUGGCCGCUCCCGUGGGACAUCGCAUAAUCUGGAGACUUAUAGGCGCAAGCAAUAGGCGGAGGUUAUGGACGGGGGAGGGAAACCGGCGGUGCAUACAAUUCCAGAUACUGAUUUCAUGCUAUUGCUCUAGCAGCCAAGUUAACAAAUACGCCUGCGCCCCGAACGAUAAUGUUUCUGCAGUGGGGACCCAGAGGUGGCUGCUUCGUGCUUGGAGCUCCCGUUUUAAUCGGGCGCAUGACCAGCACCCCACGUUGUUUUCUUUAUCCACCGCCCUAGUGGAGUGGAUGUACGCCGCCAAUAUGGCCCAAGGAAGAGACGAGCCAUCUUACCAAGAUAGCGUCGAACUUACCGACCCCUGUGUCGCUGUCCCGAGGAACGGACAAGGGAAUUACACAACAGUCCAGCUGUGGGGUUGCAAUCAUAAACGAACUCGCACAUAUUUUUCACGUGAUUGGCUAAUUCCGCCAGCAAGACGGAUAUCACAUCUAUUCACUGUAUUAUGCAAAUUUCAUCCUGUAAUCCAGUUAGGUCGAGCCACAAAUGGAAGGUUACGGACGCGUACGCCAAACUCUUUUCGUGUGGUAGGGCAGAGUCAACAUAUGGACAUCGAUCAAUUACCAUACGAAGCAGUGACUAAGCAGGUCGAGGUCAUCGUGCUAGCCGACCAAAAGCUUAGUUCGCGUGGGCAUCCAACGGACUUGAGAUUGCCGAAAGUUAGACGGACACUGAUCAGGAUCCUAAGUGGUUCUAUUGCUAUAUUAUCUAAAGUCAAUAGGCACAUGGCGAGCGUGUCGGAACCAUGGGUGGACUUGGUGACGUUGCCGGUGUUUUACGUGCUGGUGACUAAUUCGACUGACGUACAACUCAGGGAUCCUUUCGGAAGAGAAGGCUCCAGCAGGCUGGAAGUGCGGGAUCAGAGAAGUCCUUGCACAAAAUAUUUGGUUGAGGCAUUAUCUUUGGCAGGAUCAAGAUUAUACAAAACUCCCCUUACGGUCGUACAGACUGUAUAUACUAUAGGAGUUCCAUUUCGUGUCCUUAUUUAUAACACCAGGCUGUUGCGAAGCGAACGUCUGCCGUCUGAGCCACCGUCUAGGGGGAAACCUCUCUUAGAAAGCGUGCUUGAACCUCCCAGGGCUACCGACUCGGAGAUCCAGGCGGCGACGGUUCUAAAACGACUUCGAGGCAGUUGUCCCUUAAAUGUAGCGCAUCAAAUAGACUGUUUAUCCCCCCCACGACAUGACCGCAAGAAGGGGUCCCACGUCGGGACGUUUAUCCCCGUGCGAAUGGCUACUUGCGACAGCCACAUGAACGGCAAGAUGCUAGGUUGGUCCGCGGACCUGACUUUAAGCAACAAUGGUGCCCCAUUUGCAAGGCUAGCCACAGGGAAGGCCCCGCUUUGUCGUACUCUCUUGCCAAUCGUUAAUGCUUGGUUAAGGAUCAAGGUUCUCGGAUGGACCGCAAGAUGUGCGGUGACUAGAUCUGGACCUCACGCCCGCCGACUAGCUCUCCUCAGUAAUUCAAGUUCAAAAAACAGAGUGACUACUUUAGGCGACGUGCGCAUCCUCAUGUUCGCCGGAUCUAACGCUAAAGGGUUAUGUGUUCCGCUGGAACUGAGGGUAUUUGCCCCUCUUGCCCGCCCUUACUAUAGUAACCAUGCUGGUGGCUUCCGCUUAGUAUACGCCCGCUUUACAACUAGAGCUCGUUGUGCCAUGGUAUGCGUUGGUCAGAACGACUACCUAAAAGGGAAUCGGGGUUCGGCAAUUCACAGGGAGUCCAAUCUGCUCUCUUUAUACGCUCGUAGGAAUAAAGGAUCUAAACCAACAGCGCAAUACCACACCGCGGCGCUCAAAAUUGCUUCAGUGACAAAUCUAAGGUACAUUGUAUGGCCAAGACCUAAGGAUCGCCUUAAAAAACACCGUAGCGCUCGGUGGUUAUGCAGCGCCCCAAAACUUUUACUUAUCUACGGUGCAAGAGCCCUAGGAGGAGAAGUCACAUACGAACUACUUGGGAGUGUUCAGAGUGGAACGCUUCUUCCACGUACUCUCCCUCCUCGGACUUGCGGUCUACCAUUAAACACCAUGUCAGUGUAUAGUAAUAACACGUCUAGGUUGAUGGUUCCGCAUAAAUGUCCGACUCGCCGAAUUGCGUUGUGGUGGAUUAAGAAGGCGCCACCUGGAGAUUCAAUACGAGUGUGGGUAGUAAAGUCAUCGGUCAGACCUGAUUGCAAUCAGCGUCGAACCUUAGGCGAACAAAAGACGUUACCAGGGAGGCCGGAGACCUCAAGACUGCACGUCUACAGCGAGGUGACGGUAGAUGGUGCUGGCUCAACCAAGCAACUUGAUGUGCCUACACUAAAUCUAACCAACGAGUUCUGCUUGCCGGCUUCAUCCAUAGCACCGCUGAUAGUAUCCUCCCACCCAUCUGGUGUUCAACACAAUGCGUUUUGCGCGACAAUGCCCAUUCUCGACGGGUGGUUCAAGCUGAGUCAGACGUGUAAUGUUUCUAUCAGUAACGGGCACGUAGAUAGAGACGUCCCACUUAUCUACUGCUGUGCUUUUAGUAUGCCGCGCCGAACGUACGUCCAGGUCCAUCUCUGUUCGAGCCGGGGGGCUAUGACAGCGGAUCGCAACAAUCCGCACAAAGUGGGGAGAUACCCCGGCUUGUACUCUCAUGAACAGUUCACUCAUGCAUCAGCGAGCUUGCGGUCGGGCACGGAAGGUAGAAGGCCAUGCGGCUCGAAGUACAGCUACGGUUCGCGAGGGGGUCGAGAUCGUGCCAGAAGGAUGCAAGGUGGGGAAAGCGACUUAAGCUCGAGGAUCCUCAGGCUGAUCCUAUUAAUGGUCAUGCCGUACGUUACUCGUGAACGCAGUCAUUCAGGCGAGCCAUUUGGCUUAUCCCGACGGCGGCCGCUGGCAUCUAGGAUAGAGUACACUCUGCCGUUAUCCCUUGGCGGAGGUCUGGACGAGAACCAAGAGUAUCGGCAGAUCCCUGCAAGUAGGAGAUGUAGGGCUGGACUUGACAGAAUGUUGACGAGGAUAACAAUCCCACUUAGGUUGGAGGACACUCUUCUUCGUAUAGCGGGGACCGGAUCCCGAUGCGGCUUAAAAGAUUUGAUUGGCUCAUCAAGGGUGAUGUCUGCCCACCAAAUCAACAUCACGCCGCAGUCAGCAGAACGGACGUCCUUAAGUUCGGGUCUGUACUACCCCUACCAACACAAAUUUGACACAGACAGACAGCUCCGGAAGUUACCGCUUGGUGAUCUCUCUAAACUAUGUCUGGCAGUGGGAACACGUUAUAACGAGUAUCGUGUUUUGGUUCAAUGCCAUCUACUAGUUAACACUCACGGCCCGGGUUGGGAGUUCUAUAACACGAGGCAGAAUUCUGCUUUUGCUGUACAACAUAGUUUUAUAGGAAUGCGAUGCGAAUUGCUGGAAGAGAGAAUGGGCUGGAUGGCGAUAAGGGAAAACAAGUUCUCACGUCCUGACUGUACGGCAAUGCGUUGUGGAGUAAGCCGGGCGUGUAUGUUAAGCAUCGGACCUGGGCCAGAGGGGUUUGAGUUGGACCAGAAAAAGCUGCCUCACCUGCUAGAAAUGGUGCAAGUGUCACUCCGGAUUAUUCAAAUAAUGACCUGGGACGACGCCCCCCAGUUUACGGGGCGGGCCCAUCCCCUACAAGUCUGCUCGCUACACAAAGUCCAUUUCUCUACCAGCCGCUGCCGAAAGAUGUCACGAGAUAUUGUGCCGCUGCCAGCUGGGCUUGUGUCCUUGCCGGGCAGAGAGCACCCUGCAAUGGUGAUUCGUGGCUGCCUUCGCCCUACAAAGUCAGGACACCGCAGUGGAUUGAGGUGCUUCUCGAGGCCUUCGACGGUGAUGCAAAAGGACAAGGCACUAGAUGAACUCCUGAUGGCACUGUACCAGUACCACCCCAUCUCCACAUGA"

val motifPair = Pair(text, subs)

fun rosalindCountNucleotides() {
    println()
    print("Rosalind Count Nucleotides: ")
    println(rbt.countBases(dna))
}

fun rosalindDNACompliment() {
    println()
    print("Rosalind Reverse Compliment: ")
    println(rbt.reverseCompliment(dna))
}

fun rosalindGCContent() {
    val strands = rbt.stripFASTA("../data/testFileGC-Content")
    println()
    print("Rosalind GC Content: ")
    println(rbt.gcContent(strands))
}

fun rosalindHammingDistance() {
    val pairs = rbt.readHammingPair("../data/testFileHamming")
    println()
    print("Rosalind Hamming Distance: ")
    println(rbt.calculateHammingDistance(pairs))
}

fun rosalindTranscribeDNAToRNA() {
    println()
    print("Rosalind Transcribe DNA to RNA: ")
    println(rbt.transDnaToRna(dna))
}

fun rosalindMotifLocations() {
    println()
    print("Rosalind Motif Locations: ")
    println(rbt.motifLocations(motifPair))
}

fun rosalindSharedMotifs() {
    println()
    print("Rosalind Shared Motifs: ")
    val foo = rbt.stripFASTA("../data/testFileSharedMotif")
    val bar = rbt.uniqMotifs(foo)
    println(rbt.largestCommonMotif(bar))
}

fun rosalindPermutations(num: Int) {
    println()
    println("Rosalind Permutations: ")
    println((1..num).toList().permutationsNumber())
    for (x in (1..num).toList().permutations()) {
        for (perm in x) {
            print(perm.toString() + " ")
        }
        println()
    }
}

fun rosalindTranscribeRNAToProtein() {
    println()
    print("Rosalind Transcribe RNA To Protein: ")
    println(rbt.transRna2Protein(rna2protein))

}

fun rosalindProteinMotifs(inputUrl: String) {
    println()
    print("Rosalind Protein Motif Locations: ")
    for (x in rbt.proteinMotifLocations(inputUrl)) {
        print("${x} ")
    }
}

rosalindCountNucleotides()
rosalindDNACompliment()
rosalindGCContent()
rosalindHammingDistance()
rosalindTranscribeDNAToRNA()
rosalindMotifLocations()
//rosalindSharedMotifs()
rosalindPermutations(5)
rosalindTranscribeRNAToProtein()
rosalindProteinMotifs("https://www.uniprot.org/uniprot/P07204_TRBM_HUMAN.fasta")

