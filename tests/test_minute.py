from minute import read_tsv


def test_read_tsv(tmp_path):
    tsv = tmp_path / "file.tsv"
    tsv.write_text("a\t  b c\t")
    assert list(read_tsv(tsv, 3)) == [
        ["a", "b", "c"],
    ]
