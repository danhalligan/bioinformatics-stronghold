import rosalind.cli as cli
import pytest
import yaml
from unittest.mock import patch

# Dictionary of saved uniprot output
uniprot = yaml.load(open("tests/uniprot_output.yaml"), Loader=yaml.FullLoader)

# Gather registered commands
cmds = [x.name for x in cli.app.registered_commands]


# Run each command with a test file and snapshot output
@pytest.mark.parametrize("fun", cmds)
def test_cli_function(capfd, snapshot, fun):
    with patch("rosalind.rosalind.uniprot_output", side_effect=uniprot.get):
        getattr(cli, fun)(f"tests/data/test_{fun}.txt")
    out, err = capfd.readouterr()
    snapshot.assert_match(out)


def test_gaff(capfd, snapshot):
    for i in range(1, 4):
        out = cli.gaff(f"tests/data/test_gaff{i}.txt")
        out, err = capfd.readouterr()
        snapshot.assert_match(out)
