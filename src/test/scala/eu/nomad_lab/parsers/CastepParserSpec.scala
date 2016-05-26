package eu.nomad_lab.parsers

import org.specs2.mutable.Specification

object CastepParserSpec extends Specification {
  "CastepParserTest" >> {
    "test with json-events" >> {
      ParserRun.parse(CastepParser, "parsers/castep/test/examples/Si2_bandstructure/B3LYP/Si2.castep_b_sp_v_1", "json-events") must_== ParseResult.ParseSuccess
    }
    "test with json" >> {
      ParserRun.parse(CastepParser, "parsers/castep/test/examples/Si2_bandstructure/B3LYP/Si2.castep_b_sp_v_1", "json") must_== ParseResult.ParseSuccess
    }
  }
}

object CastepParserSpec1 extends Specification {
  "CastepParserTest" >> {
    "test with json-events" >> {
      ParserRun.parse(CastepParser, "parsers/castep/test/examples/KCl-stress.castep", "json-events") must_== ParseResult.ParseSuccess
    }
    "test with json" >> {
      ParserRun.parse(CastepParser, "parsers/castep/test/examples/KCl-stress.castep", "json") must_== ParseResult.ParseSuccess
    }
  }
}

object CastepParserSpec extends Specification {
  "CastepParserTest" >> {
    "test with json-events" >> {
      ParserRun.parse(CastepParser, "parsers/castep/test/examples/Si2_singlepoint/LDA/Si2.castep_v_2", "json-events") must_== ParseResult.ParseSuccess
    }
    "test with json" >> {
      ParserRun.parse(CastepParser, "parsers/castep/test/examples/Si2_singlepoint/LDA/Si2.castep_v_2", "json") must_== ParseResult.ParseSuccess
    }
  }
}

