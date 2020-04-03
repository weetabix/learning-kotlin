//for (i in 2..num-1) {
//    if (i % 2 != 0) {
//        if (num % i == 0) {
//          //  println(i)
//            for (z in 2..i-1) {
//                if (z % i == 0) {
//                    if (i % z == 0) {
//                        println(z)
//                    }
//}
//}
//        }
//    }
//}


var num = 13195


class Problem3 {

    fun isPrime(inum: Int) {
        for (x in 1..inum) {
            if (inum % x == 0) {
                println(x)
            }
        }
    }
}

val foo = Problem3()

foo.isPrime(num)

