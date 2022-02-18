from pathlib import Path

from minute import (
    read_tsv,
    read_libraries,
    read_scaling_groups,
    Library,
    libraries_unused_in_groups,
)


def test_read_tsv(tmp_path):
    tsv = tmp_path / "file.tsv"
    tsv.write_text("a\t  b c\t")
    assert list(read_tsv(tsv, 3)) == [
        ["a", "b", "c"],
    ]


def test_libraries_unused_in_groups_tsv():
    libraries = list(read_libraries(Path("testdata") / "libraries.tsv"))
    groups = list(read_scaling_groups(Path("testdata") / "groups.tsv", libraries))

    assert len(libraries_unused_in_groups(libraries, groups)) == 0

    unused_lib = Library("unused")
    libraries.append(unused_lib)
    unused = libraries_unused_in_groups(libraries, groups)
    assert len(unused) == 1
    assert unused[0] is unused_lib
