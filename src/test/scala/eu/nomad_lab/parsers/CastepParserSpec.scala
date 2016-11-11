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

object CastepParserSpec2 extends Specification {
  "CastepParserTest" >> {
    "test with json-events" >> {
      ParserRun.parse(CastepParser, "parsers/castep/test/examples/Si2_singlepoint/LDA/Si2.castep_v_2", "json-events") must_== ParseResult.ParseSuccess
    }
    "test with json" >> {
      ParserRun.parse(CastepParser, "parsers/castep/test/examples/Si2_singlepoint/LDA/Si2.castep_v_2", "json") must_== ParseResult.ParseSuccess
    }
  }
}
object CastepParserSpec3 extends Specification {
  "CastepParserTest" >> {
    "test with json-events" >> {
      ParserRun.parse(CastepParser, "parsers/castep/test/examples/Si2_singlepoint/LDA/Spin_polarised/Si2_sp.castep", "json-events") must_== ParseResult.ParseSuccess
    }
    "test with json" >> {
      ParserRun.parse(CastepParser, "parsers/castep/test/examples/Si2_singlepoint/LDA/Spin_polarised/Si2_sp.castep", "json") must_== ParseResult.ParseSuccess
    }
  }
}
object CastepParserSpec4 extends Specification {
  "CastepParserTest" >> {
    "test with json-events" >> {
      ParserRun.parse(CastepParser, "parsers/castep/test/examples/si2-bfgs.castep", "json-events") must_== ParseResult.ParseSuccess
    }
    "test with json" >> {
      ParserRun.parse(CastepParser, "parsers/castep/test/examples/si2-bfgs.castep", "json") must_== ParseResult.ParseSuccess
    }
  }
}
object CastepParserSpec5 extends Specification {
  "CastepParserTest" >> {
    "test with json-events" >> {
      ParserRun.parse(CastepParser, "parsers/castep/test/examples/MD/Si8-md-NPT.castep", "json-events") must_== ParseResult.ParseSuccess
    }
    "test with json" >> {
      ParserRun.parse(CastepParser, "parsers/castep/test/examples/MD/Si8-md-NPT.castep", "json") must_== ParseResult.ParseSuccess
    }
  }
}
object CastepParserSpec6 extends Specification {
  "CastepParserTest" >> {
    "test with json-events" >> {
      ParserRun.parse(CastepParser, "parsers/castep/test/examples/NH3-TS-PhonGamma-400.castep", "json-events") must_== ParseResult.ParseSuccess
    }
    "test with json" >> {
      ParserRun.parse(CastepParser, "parsers/castep/test/examples/NH3-TS-PhonGamma-400.castep", "json") must_== ParseResult.ParseSuccess
    }
  }
}

object CastepParserSpec7 extends Specification {
  "CastepParserTest" >> {
    "test with json-events" >> {
      ParserRun.parse(CastepParser, "parsers/castep/test/examples/Si2-cellopt-aniso.castep", "json-events") must_== ParseResult.ParseSuccess
    }
    "test with json" >> {
      ParserRun.parse(CastepParser, "parsers/castep/test/examples/Si2-cellopt-aniso.castep", "json") must_== ParseResult.ParseSuccess
    }
  }
}

object CastepParserSpec8 extends Specification {
  "CastepParserTest" >> {
    "test with json-events" >> {
      ParserRun.parse(CastepParser, "parsers/castep/test/examples/TS_search/h2-lst.castep", "json-events") must_== ParseResult.ParseSuccess
    }
    "test with json" >> {
      ParserRun.parse(CastepParser, "parsers/castep/test/examples/TS_search/h2-lst.castep", "json") must_== ParseResult.ParseSuccess
    }
  }
}
object CastepParserSpec9 extends Specification {
  "CastepParserTest" >> {
    "test with json-events" >> {
      ParserRun.parse(CastepParser, "parsers/castep/test/examples/N2-LDA-geom.castep", "json-events") must_== ParseResult.ParseSuccess
    }
    "test with json" >> {
      ParserRun.parse(CastepParser, "parsers/castep/test/examples/N2-LDA-geom.castep", "json") must_== ParseResult.ParseSuccess
    }
  }
}
object CastepParserSpec10 extends Specification {
  "CastepParserTest" >> {
    "test with json-events" >> {
      ParserRun.parse(CastepParser, "parsers/castep/test/examples/Surf/surf.castep", "json-events") must_== ParseResult.ParseSuccess
    }
    "test with json" >> {
      ParserRun.parse(CastepParser, "parsers/castep/test/examples/Surf/surf.castep", "json") must_== ParseResult.ParseSuccess
    }
  }
}
object CastepParserSpec11 extends Specification {
  "CastepParserTest" >> {
    "test with json-events" >> {
      ParserRun.parse(CastepParser, "parsers/castep/test/examples/test.castep", "json-events") must_== ParseResult.ParseSuccess
    }
    "test with json" >> {
      ParserRun.parse(CastepParser, "parsers/castep/test/examples/test.castep", "json") must_== ParseResult.ParseSuccess
    }
  }
}
object CastepParserSpec12 extends Specification {
  "CastepParserTest" >> {
    "test with json-events" >> {
      ParserRun.parse(CastepParser, "parsers/castep/test/examples/benzene_cell/benzene_cell.castep", "json-events") must_== ParseResult.ParseSuccess
    }
    "test with json" >> {
      ParserRun.parse(CastepParser, "parsers/castep/test/examples/benzene_cell/benzene_cell.castep", "json") must_== ParseResult.ParseSuccess
    }
  }
}
object CastepParserSpec13 extends Specification {
  "CastepParserTest" >> {
    "test with json-events" >> {
      ParserRun.parse(CastepParser, "parsers/castep/test/examples/Gph-H/Gph-H.castep", "json-events") must_== ParseResult.ParseSuccess
    }
    "test with json" >> {
      ParserRun.parse(CastepParser, "parsers/castep/test/examples/Gph-H/Gph-H.castep", "json") must_== ParseResult.ParseSuccess
    }
  }
}
object CastepParserSpec14 extends Specification {
  "CastepParserTest" >> {
    "test with json-events" >> {
      ParserRun.parse(CastepParser, "parsers/castep/test/examples/CHchair/CH-chair_bands_DOS.castep", "json-events") must_== ParseResult.ParseSuccess
    }
    "test with json" >> {
      ParserRun.parse(CastepParser, "parsers/castep/test/examples/CHchair/CH-chair_bands_DOS.castep", "json") must_== ParseResult.ParseSuccess
    }
  }
}
object CastepParserSpec15 extends Specification {
  "CastepParserTest" >> {
    "test with json-events" >> {
      ParserRun.parse(CastepParser, "parsers/castep/test/examples/test2/CF-boat.castep", "json-events") must_== ParseResult.ParseSuccess
    }
    "test with json" >> {
      ParserRun.parse(CastepParser, "parsers/castep/test/examples/test2/CF-boat.castep", "json") must_== ParseResult.ParseSuccess
    }
  }
}
object CastepParserSpec16 extends Specification {
  "CastepParserTest" >> {
    "test with json-events" >> {
      ParserRun.parse(CastepParser, "parsers/castep/test/examples/c4h3/C4H3.castep", "json-events") must_== ParseResult.ParseSuccess
    }
    "test with json" >> {
      ParserRun.parse(CastepParser, "parsers/castep/test/examples/c4h3/C4H3.castep", "json") must_== ParseResult.ParseSuccess
    }
  }
}
object CastepParserSpec17 extends Specification {
  "CastepParserTest" >> {
    "test with json-events" >> {
      ParserRun.parse(CastepParser, "parsers/castep/test/examples/CH-boatbands.castep", "json-events") must_== ParseResult.ParseSuccess
    }
    "test with json" >> {
      ParserRun.parse(CastepParser, "parsers/castep/test/examples/CH-boatbands.castep", "json") must_== ParseResult.ParseSuccess
    }
  }
}
object CastepParserSpec18 extends Specification {
  "CastepParserTest" >> {
    "test with json-events" >> {
      ParserRun.parse(CastepParser, "parsers/castep/test/examples/ch-chair2/CH-chairc04x00.castep", "json-events") must_== ParseResult.ParseSuccess
    }
    "test with json" >> {
      ParserRun.parse(CastepParser, "parsers/castep/test/examples/ch-chair2/CH-chairc04x00.castep", "json") must_== ParseResult.ParseSuccess
    }
  }
}

object CastepParserSpec19 extends Specification {
  "CastepParserTest" >> {
    "test with json-events" >> {
      ParserRun.parse(CastepParser, "parsers/castep/test/examples/Ag_2x4/surf.castep", "json-events") must_== ParseResult.ParseSuccess
    }
    "test with json" >> {
      ParserRun.parse(CastepParser, "parsers/castep/test/examples/Ag_2x4/surf.castep", "json") must_== ParseResult.ParseSuccess
    }
  }
}

object CastepParserSpec20 extends Specification {
  "CastepParserTest" >> {
    "test with json-events" >> {
      ParserRun.parse(CastepParser, "parsers/castep/test/examples/5Sb2Te3only426Sym/5Sb2Te3only426Sym.castep", "json-events") must_== ParseResult.ParseSuccess
    }
    "test with json" >> {
      ParserRun.parse(CastepParser, "parsers/castep/test/examples/5Sb2Te3only426Sym/5Sb2Te3only426Sym.castep", "json") must_== ParseResult.ParseSuccess
    }
  }
}

