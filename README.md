# `rosalind`

Solutions to the Rosalind problems in python as a command line app. My solutions
may be brute force :relaxed:

<http://rosalind.info/problems/>

I'm trying not to use Biopython (whilst also not attempting to rewrite this from
scratch!)

Some build details that I'm playing with while putting this together:

- CLI built using [Poetry] & [Typer] based on this [blog post][pluralsight]
- [GitLab CI] based on this medium [blog post][medium]
- [black] for code formatting
- [flake8] for code linting
- Pre-commit hook for formatting with black using [pre-commit] based on
  [black docs]
- Coverage and coverage reports with [pytest-cov]
- Snapshots of test output with [snapshottest]

## Install

You can install with pip:

```bash
git clone https://gitlab.fiosgenomics.com/dan/rosalind rosalind
pip3 install rosalind
```

## Examples

``` bash
rosalind dna rosalind_dna.txt
```

## Testing

You can run tests using poetry. For example, to test a function with the test
data, run something like:

``` bash
poetry run rosalind fib tests/data/test_fib.txt
```

To run all tests and calculate coverage run:

``` bash
poetry run pytest --cov rosalind
```

[Poetry]: https://python-poetry.org/
[Typer]: https://typer.tiangolo.com/
[pluralsight]: https://www.pluralsight.com/tech-blog/python-cli-utilities-with-poetry-and-typer/
[black]: https://black.readthedocs.io/en/stable/index.html
[flake8]: https://flake8.pycqa.org/en/latest/
[GitLab CI]: https://docs.gitlab.com/ee/ci/
[medium]: https://medium.com/@paweldudzinski/python-applications-continuous-integration-with-poetry-and-gitlab-pipelines-ac539888251a
[pre-commit]: https://pre-commit.com/
[black docs]: https://black.readthedocs.io/en/stable/version_control_integration.html
[pytest-cov]: https://pypi.org/project/pytest-cov/
[snapshottest]: https://pypi.org/project/snapshottest/
