package eu.nomad_lab.parsers
import eu.nomad_lab.DefaultPythonInterpreter
import org.{json4s => jn}

object CastepParser extends SimpleExternalParserGenerator(
  name = "CastepParser",
  parserInfo = jn.JObject(
    ("name" -> jn.JString("CastepParser")) ::
      ("version" -> jn.JString("1.0")) :: Nil),
  mainFileTypes = Seq("text/.*"),
  mainFileRe = """\s\|\s*CCC\s*AA\s*SSS\s*TTTTT\s*EEEEE\s*PPPP\s*\|\s*""".r,
  cmd = Seq(DefaultPythonInterpreter.python2Exe(), "${envDir}/parsers/castep/parser/parser-castep/CastepParser.py",
    "--uri", "${mainFileUri}", "${mainFilePath}"),
  resList = Seq(
    "parser-castep/CastepParser.py",
    "parser-castep/CastepBandParser.py",
    "parser-castep/CastepCellParser.py",
    "parser-castep/CastepCommon.py",
    "parser-castep/setup_paths.py",
    "nomad_meta_info/common.nomadmetainfo.json",
    "nomad_meta_info/meta_types.nomadmetainfo.json",
    "nomad_meta_info/castep.nomadmetainfo.json"
  ) ++ DefaultPythonInterpreter.commonFiles(),
  dirMap = Map(
    "parser-castep" -> "parsers/castep/parser/parser-castep",
    "nomad_meta_info" -> "nomad-meta-info/meta_info/nomad_meta_info"
  ) ++ DefaultPythonInterpreter.commonDirMapping()
)
