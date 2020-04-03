class Hamming() {

    var hamm_count = 0

    fun calculateHammingDistance(fileName: String) {


        for (b in 0..stringA.length - 1) {
            if (stringA[b] != stringB[b])
                hamm_count++
        }
        println(hamm_count)
    }
}

val foo = Hamming()
foo.calculateHammingDistance("../testFileHamming")