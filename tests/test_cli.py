from pathlib import Path

from atomigraph.cli import main

EXAMPLES = Path(__file__).resolve().parents[1] / "examples"


def test_cli_finds_and_plots_reactions(tmp_path, monkeypatch, capsys):
    monkeypatch.chdir(tmp_path)
    infile = str(EXAMPLES / "01_PE_chain_basic" / "pe_chain.*.data")

    rc = main(["-i", infile, "-f", "lammps_data", "-b", "plots", "--plot-format", "png"])

    assert rc == 0
    assert "5 reaction(s) found" in capsys.readouterr().out
    assert len(list((tmp_path / "plots").iterdir())) == 5


def test_cli_accepts_multiple_files(tmp_path, monkeypatch, capsys):
    # unquoted glob: the shell passes each file separately
    monkeypatch.chdir(tmp_path)
    files = sorted(str(p) for p in (EXAMPLES / "03_PEEK_multiple_files").glob("*.dump"))

    rc = main(["-i", *files, "-a", "1:C,2:H,3:H,4:O,5:O,6:O,7:O,8:O", "--no-plot"])

    assert rc == 0
    assert "1 reaction(s) found" in capsys.readouterr().out


def test_cli_missing_input_file(tmp_path, monkeypatch, capsys):
    monkeypatch.chdir(tmp_path)

    rc = main(["-i", "does_not_exist.dump"])

    assert rc == 1
    assert "does_not_exist.dump" in capsys.readouterr().err
