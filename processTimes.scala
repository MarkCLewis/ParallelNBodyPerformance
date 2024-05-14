import java.io.File

val lineRegex = """(\d+) (\d+)|real\s+(\d+)m(\d+.\d+)s""".r

def stats(times: List[Double]): (Double, Double) = {
  val avg = times.sum / times.length
  val std = math.sqrt(times.map(t => (t - avg)*(t - avg)).sum / times.length)
  (avg, std)
}

@main def processTimes(topDir: String): Unit = {
  val dirs = new File(topDir).listFiles().filter(_.isDirectory()).sortBy(_.getName())
  var id = 0
  for (dir <- dirs) {
    val files = dir.list().filter(name => name.startsWith("time") && name.endsWith(".txt"))
    for (file <- files) {
      println(s"$dir/$file")
      var threads = 0
      var bodies = 0
      var secs: List[Double] = Nil
      val source = io.Source.fromFile(s"$dir/$file")
      val lines = source.getLines()
      for (lineRegex(n, t, m, s) <- lines) {
        if (n != null) {
          if (secs.nonEmpty) {
            val (avg, std) = stats(secs)
            println(s"$id $bodies $threads $avg $std")
            secs = Nil
          }
          bodies = n.toInt
          threads = t.toInt
          if (threads % 2 == 1) threads += 1
        } else {
          secs ::= 60*m.toDouble + s.toDouble
        }
      }
      val (avg, std) = stats(secs)
      println(s"$id $bodies $threads $avg $std")
      source.close()
      id += 1
    }
  }
}
