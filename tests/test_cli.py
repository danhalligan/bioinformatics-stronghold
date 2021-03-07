import rosalind.cli as cli
import pytest


# Gather registered commands
cmds = [x.name for x in cli.app.registered_commands]


# Run each command with a test file and snapshot output
@pytest.mark.parametrize("fun", cmds)
def test_cli_function(capfd, snapshot, fun):
    getattr(cli, fun)(f"tests/data/test_{fun}.txt")
    out, err = capfd.readouterr()
    snapshot.assert_match(out)
