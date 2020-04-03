class Euler2 {
    var total: Int = 0
    fun fibonacciGenerator(limit: Int = 1000000): Int {
        var x: Int = 0
        var y: Int = 1
        var z: Int
        while (x < limit) {

            total = total + evenFilter(x)
            z = y
            y = x + y
            x = z
        }
        return total
    }

    fun evenFilter(num: Int): Int {
        return if (num % 2 == 0) num else 0
    }
}

var res = Euler2()
var fib = res.fibonacciGenerator(4000000)
println(fib)