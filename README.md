# Rosalind solver

Solutions to the Rosalind problems in python as a command line app. My solutions
may be brute force :relaxed:

<http://rosalind.info/problems/>

I'm trying not to use Biopython (whilst also not attempting to rewrite this from
scratch!)

Some build details that I'm playing with while putting this together:

- CLI built using Poetry & Typer based on this blog post.
  <https://www.pluralsight.com/tech-blog/python-cli-utilities-with-poetry-and-typer/>
- [GitLab CI] based on this [medium blog post][ci-post]
- [Black] for code formatting
- Pre-commit hook for formatting with black using [pre-commit] based on 
  [black docs]

## Install

You can install with pip:

```bash
git clone https://gitlab.fiosgenomics.com/rosalind/dan rosalind-solver
pip3 install rosalind-solver
```

[Black]: https://black.readthedocs.io/en/stable/index.html
[GitLab CI]: https://docs.gitlab.com/ee/ci/
[ci-post]: https://medium.com/@paweldudzinski/python-applications-continuous-integration-with-poetry-and-gitlab-pipelines-ac539888251a
[pre-commit]: https://pre-commit.com/
[black docs]: https://black.readthedocs.io/en/stable/version_control_integration.html
