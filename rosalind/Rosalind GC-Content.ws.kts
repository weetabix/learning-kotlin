import java.io.BufferedReader
import java.io.File

class Problem4() {
    var gc_result = Pair("None", 0.0)
    fun gcContent(fileName: String) {

        val bufferedReader: BufferedReader = File(fileName).bufferedReader()
        val inputString = bufferedReader.use { it.readText() }
        val samples = inputString.split(">").toTypedArray()
        for (i in samples.indices) {
            var gc_count = 0
            val rawsample = samples[i].split("\n").toTypedArray()
            val name = rawsample[0]
            val strand = samples[i].split("\n").toTypedArray().drop(1).joinToString(separator = "")
            for (base in 0..strand.length - 1) {
                when (strand[base].toString()) {
                    "C", "G" -> gc_count++
                }
            }
            var gc_perc = gc_count.toDouble() / strand.length * 100
            if (gc_perc > gc_result.second) {
                gc_result = Pair(name, gc_perc)
            }
        }
        println(gc_result)
    }
}

val foo = Problem4()
foo.gcContent("../testFileGC-Content")
