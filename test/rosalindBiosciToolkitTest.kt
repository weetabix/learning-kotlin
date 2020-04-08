import org.junit.jupiter.api.Test

import org.junit.jupiter.api.Assertions.*

internal class rosalindBiosciToolkitTest {

    @Test
    fun file2string() {
    }

    @Test
    fun countBases() {
        val ans = rosalindBiosciToolkit().countBases("ACGT")
        assertEquals(ans["A"], 1)
        assertEquals(ans["C"], 1)
        assertEquals(ans["G"], 1)
        assertEquals(ans["T"], 1)
    }

    @Test
    fun reverseCompliment() {
    }

    @Test
    fun stripFASTA() {
    }

    @Test
    fun readHammingPair() {
    }

    @Test
    fun gcContent() {
    }

    @Test
    fun calculateHammingDistance() {
    }

    @Test
    fun transDnaToRna() {
    }

    @Test
    fun motifLocations() {
    }

    @Test
    fun uniqMotifs() {
    }

    @Test
    fun largestCommonMotif() {
    }

    @Test
    fun coLargestCommonMotif() {
    }

    @Test
    fun transRna2Protein() {
    }
}