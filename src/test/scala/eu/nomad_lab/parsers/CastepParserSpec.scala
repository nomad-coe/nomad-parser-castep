package eu.nomad_lab.parsers

import org.specs2.mutable.Specification



object CastepParserSpec extends Specification {
  "CastepParserTest" >> {
    "test with json-events" >> {
      ParserRun.parse(CastepParser, "parsers/castep/test/examples/Si2.castep_v_1", "json-events") must_== ParseResult.ParseSuccess
    }
    "test with json" >> {
      ParserRun.parse(CastepParser, "parsers/castep/test/examples/Si2.castep_v_1", "json") must_== ParseResult.ParseSuccess
    }
  }
}
