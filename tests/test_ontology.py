from __future__ import annotations

from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

from cafa.ontology import (
    ancestors_of,
    canonicalize_go_id,
    filter_terms_to_subontology,
    parents_of,
    propagate_scores,
    propagate_terms,
    read_go_obo,
    subontology_terms,
    terms_of_interest,
)


TEST_OBO = """format-version: 1.2
ontology: go
name: mini-go
data-version: releases/2025-06-01

[Term]
id: GO:0003674
name: molecular_function
namespace: molecular_function

[Term]
id: GO:0008150
name: biological_process
namespace: biological_process

[Term]
id: GO:0005575
name: cellular_component
namespace: cellular_component

[Term]
id: GO:0000001
name: parent term
namespace: molecular_function
def: "parent"
is_a: GO:0003674 ! molecular_function

[Term]
id: GO:0000002
name: child term
namespace: molecular_function
alt_id: GO:0009002
def: "child"
is_a: GO:0000001 ! parent term
relationship: part_of GO:0003674 ! molecular_function

[Term]
id: GO:0000003
name: obsolete term
namespace: molecular_function
is_obsolete: true
is_a: GO:0000001 ! parent term

[Term]
id: GO:0000100
name: process term
namespace: biological_process
is_a: GO:0008150 ! biological_process
"""


class OntologyTests(unittest.TestCase):
    def setUp(self) -> None:
        self.tmpdir = TemporaryDirectory()
        self.obo_path = Path(self.tmpdir.name) / "mini.obo"
        self.obo_path.write_text(TEST_OBO, encoding="utf-8")
        self.ontology = read_go_obo(self.obo_path)

    def tearDown(self) -> None:
        self.tmpdir.cleanup()

    def test_read_go_obo_extracts_release_and_ignores_obsolete_terms(self) -> None:
        self.assertEqual(self.ontology.release, "2025-06-01")
        self.assertIn("GO:0000001", self.ontology.terms)
        self.assertNotIn("GO:0000003", self.ontology.terms)

    def test_canonicalize_go_id_uses_alt_id(self) -> None:
        self.assertEqual(canonicalize_go_id(self.ontology, "GO:0009002"), "GO:0000002")
        self.assertEqual(canonicalize_go_id(self.ontology, "GO:9999999"), "GO:9999999")

    def test_parents_and_ancestors(self) -> None:
        self.assertEqual(parents_of(self.ontology, "GO:0000002"), ("GO:0000001", "GO:0003674"))
        self.assertEqual(ancestors_of(self.ontology, "GO:0000002"), frozenset({"GO:0000001", "GO:0003674"}))

    def test_propagate_terms(self) -> None:
        propagated = propagate_terms(self.ontology, {"GO:0009002"})
        self.assertEqual(propagated, frozenset({"GO:0000002", "GO:0000001", "GO:0003674"}))

    def test_propagate_scores(self) -> None:
        propagated = propagate_scores(self.ontology, {"GO:0000002": 0.4, "GO:0000001": 0.2})
        self.assertEqual(propagated["GO:0000002"], 0.4)
        self.assertEqual(propagated["GO:0000001"], 0.4)
        self.assertEqual(propagated["GO:0003674"], 0.4)

    def test_subontology_filtering_and_terms_of_interest(self) -> None:
        self.assertEqual(
            subontology_terms(self.ontology, "MF"),
            frozenset({"GO:0003674", "GO:0000001", "GO:0000002"}),
        )
        self.assertEqual(
            terms_of_interest(self.ontology, "MF"),
            ("GO:0000001", "GO:0000002"),
        )
        self.assertEqual(
            filter_terms_to_subontology(self.ontology, {"GO:0000100", "GO:0009002"}, "MF"),
            frozenset({"GO:0000002"}),
        )


if __name__ == "__main__":
    unittest.main()
