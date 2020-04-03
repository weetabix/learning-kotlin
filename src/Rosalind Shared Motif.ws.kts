import java.io.BufferedReader
import java.io.File

class sharedMotif() {

    fun stripFASTA(fileName: String): List<Pair<String, String>> {

        var strands = listOf<Pair<String, String>>()

        val bufferedReader: BufferedReader = File(fileName).bufferedReader()
        val inputString = bufferedReader.use { it.readText() }
        val samples = inputString.split(">").toTypedArray()
        for (i in samples.indices) {
            strands += (Pair(
                samples[i].split("\n").toTypedArray()[0],
                samples[i].split("\n").toTypedArray().drop(1).joinToString(separator = "")
            ))
        }

        return strands.filterNot { it == Pair("", "") }
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
}

//val foo = sharedMotif()
//var derp = foo.stripFASTA("../testFileSharedMotif")
//var herp = foo.uniqMotifs(derp)
//var lerp = foo.largestCommonMotif(herp)
//for (x in lerp){
//  print("$x ")
//}






