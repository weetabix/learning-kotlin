class Euler1 {
    var total: Int = 0

    fun determineTotal() {
        for (x in 0..999) {
            if (x % 3 == 0 || x % 5 == 0) {
                total = total + x
            }
        }
        println(total)
    }
}

var prob = Problem1()
prob.determineTotal()
