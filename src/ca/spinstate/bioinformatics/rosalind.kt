package ca.spinstate.bioinformatics

import java.io.BufferedReader
import java.io.File
import java.net.URL

class Rosalind() {

    // Resources for the Rosalind Code Exersizes

    internal fun uniprotGetString(protID: String): String {
        val outString = URL("https", "www.uniprot.org", "/uniprot/$protID.fasta").readText()
        return outString
    }

    internal fun fileGetString(fileName: String): String {
        val bufferedReader: BufferedReader = File(fileName).bufferedReader()
        val outString = bufferedReader.use { it.readText() }
        return outString
    }

    internal fun seqMotifToRegex(inputMotif: String): Regex {
        val foo = inputMotif.replace("\\{".toRegex(), "[^")
        val bar = foo.replace("\\}".toRegex(), "]")
        val expr = "${bar.substring(0, 1)}(?=${bar.removeRange(0, 1)})"
        print("$expr ")
        return expr.toRegex()
    }

    fun stripFASTA(inString: String): List<Pair<String, String>> {
        var strands = listOf<Pair<String, String>>()
        val samples = inString.split(">").toTypedArray()
        for (i in samples.indices) {
            strands += (Pair(
                samples[i].split("\n").toTypedArray()[0],
                samples[i].split("\n").toTypedArray().drop(1).joinToString(separator = "")
            ))
        }
        return strands.filterNot { it == Pair("", "") }
    }

    fun readHammingPair(fileName: String): Pair<String, String> {
        val bufferedReader: BufferedReader = File(fileName).bufferedReader()
        val inputString =
            bufferedReader.use { it.readText() } // TODO <-----------------------------------------------------------------------Here
        val strings = inputString.split("\n").toTypedArray()
        val outPair = Pair(strings[0], strings[1])
        return outPair
    }

    fun countBases(strand: String): Map<String, Int> {
        val basesMap = mutableMapOf("A" to 0, "C" to 0, "G" to 0, "T" to 0)
        for (base in strand.chunked(1)) {
            when (base) {
                "A" -> basesMap.computeIfPresent("A") { _, v -> v + 1 }
                "G" -> basesMap.computeIfPresent("G") { _, v -> v + 1 }
                "C" -> basesMap.computeIfPresent("C") { _, v -> v + 1 }
                "T" -> basesMap.computeIfPresent("T") { _, v -> v + 1 }
                else -> println("Glitch " + base)
            }
        }
        return basesMap
    }

    fun reverseCompliment(strand: String): String {
        val r_strand = strand.reversed()
        var result = ""
        for (base in 0..r_strand.length - 1) {
            when (r_strand[base].toString()) {
                "T" -> result += "A"
                "A" -> result += "T"
                "C" -> result += "G"
                "G" -> result += "C"
                else -> print("Glitch: " + r_strand[base])
            }
        }
        return result
    }

    fun gcContent(strandList: List<Pair<String, String>>): Pair<String, Double> {
        var gc_result = Pair("None", 0.0)
        for (strand in strandList) {
            var gc_count = 0

            for (base in 0..strand.second.length - 1) {
                when (strand.second[base].toString()) {
                    "C", "G" -> gc_count++
                }
            }
            val gc_perc = gc_count.toDouble() / strand.second.length * 100
            if (gc_perc > gc_result.second) {
                gc_result = Pair(strand.first, gc_perc)
            }
        }
        return gc_result
    }

    fun calculateHammingDistance(stringPair: Pair<String, String>): Int {
        var hammCount = 0
        for (b in 0..stringPair.first.length - 1) {
            if (stringPair.first[b] != stringPair.second[b])
                hammCount++
        }
        return hammCount
    }

    fun transDnaToRna(strand: String): String {
        var result = ""
        for (base in 0..strand.length - 1) {
            when (strand[base].toString()) {
                "T" -> result += "U"
                else -> result += strand[base]
            }
        }
        return result
    }

    fun motifLocations(motifPair: Pair<String, String>): List<Int> {
        // Pair(strand, subs)
        var locList: List<Int> = mutableListOf()
        val strand = motifPair.first
        val subs = motifPair.second
        val expr = "${subs.substring(0, 1)}(?=${subs.removeRange(0, 1)})"
        print("$expr ")
        val regex = expr.toRegex()
        val results = regex.findAll(strand)
        for (x in results) {
            locList = locList + (x.range.start + 1)
        }
        return locList
    }

    fun uniqMotifs(pairList: List<Pair<String, String>>): List<List<String>> {
        var uniqMotif = mutableListOf<List<String>>()
        val compStrand = pairList[0].second
        for (strandPair in pairList.drop(1)) {
            var temp_results = mutableListOf<String>()
            for (x in 0..strandPair.second.length - 1) {
                for (y in 0..strandPair.second.length - x) {
                    if (strandPair.second.substring(x, x + y).isNotEmpty() &&
                        compStrand.contains(strandPair.second.substring(x, x + y)) == true
                    ) {
                        temp_results.add(strandPair.second.substring(x, x + y))
                    }
                }
            }
            uniqMotif.add(temp_results.distinct())
        }
        return uniqMotif
    }

    fun largestCommonMotif(motifs: List<List<String>>): List<String> {

        var motLen: Int = 0
        var motSet: Set<String> = motifs[0].toSet()
        for (testMotif in 1..motifs.size - 1) {
            var motSet_new = motifs[testMotif].intersect(motSet)
            motSet = motSet_new
        }
        for (item in motSet) {
            if (item.length > motLen) {
                motLen = item.length
            }
        }
        return motSet.filter { it.length == motLen }
    }

    fun coLargestCommonMotif(motifs: List<List<String>>): List<String> {

        var motLen: Int = 0
        var motSet: Set<String> = motifs[0].toSet()
        for (testMotif in 1..motifs.size - 1) {

            val motSetNew = motifs[testMotif].intersect(motSet)
            motSet = motSetNew
        }
        for (item in motSet) {
            if (item.length > motLen) {
                motLen = item.length
            }
        }
        return motSet.filter { it.length == motLen }
    }

    fun transRna2Protein(strand: String): String {
        var result = ""
        for (base in strand.chunked(3)) {
            when (base) {
                "UUA", "CUA", "UUG", "CUG", "CUU", "CUC" -> result += "L"
                "CGU", "CGC", "CGA", "AGA", "CGG", "AGG" -> result += "R"
                "UCU", "UCC", "UCA", "UCG", "AGU", "AGC" -> result += "S"
                "GUA", "GUU", "GUC", "GUG" -> result += "V"
                "ACG", "ACA", "ACU", "ACC" -> result += "T"
                "CCG", "CCA", "CCC", "CCU" -> result += "P"
                "GCU", "GCC", "GCA", "GCG" -> result += "A"
                "GGG", "GGA", "GGC", "GGU" -> result += "G"
                "AUC", "AUU", "AUA" -> result += "I"
                "UAA", "UAG", "UGA" -> result += "Stop"
                "UUC", "UUU" -> result += "F"
                "UAU", "UAC" -> result += "Y"
                "CAC", "CAU" -> result += "H"
                "AAC", "AAU" -> result += "N"
                "GAU", "GAC" -> result += "D"
                "AAA", "AAG" -> result += "K"
                "GAG", "GAA" -> result += "E"
                "CAG", "CAA" -> result += "Q"
                "UGC", "UGU" -> result += "C"
                "AUG" -> result += "M"
                "UGG" -> result += "W"
                else -> result += "---"
            }
        }
        return result
    }

    fun proteinMotifLocations(inUrlIDs: String): List<Pair<String, List<Int>>> {                 //List<Int> {
        var pMotifLoc: MutableList<Pair<String, List<Int>>> = mutableListOf()
        val regex = "N(?=([^P][ST][^P]))".toRegex()
        for (url: String in inUrlIDs.split(" ")) {
            val strands = stripFASTA(uniprotGetString(url))
            for (strand in strands) {
                var locList: MutableList<Int> = mutableListOf<Int>()
                val results = regex.findAll(strand.second)
                for (x in results) {
                        locList.add(x.range.start + 1)
                }
                if (locList.isNotEmpty()) {
                    pMotifLoc.add(Pair(url, locList))
                }
            }
        }
        return pMotifLoc
    }
}
