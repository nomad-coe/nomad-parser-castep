package eu.nomad_lab.parsers

import org.specs2.mutable.Specification



object CastepParserSpec extends Specification {
  "CastepParserTest" >> {
    examplesBlock {
      ParserRun.parse(CastepParser,"parsers/castep/test/examples/Si2.castep","json-events") must_== ParseResult.ParseSuccess
    }
  }
}
