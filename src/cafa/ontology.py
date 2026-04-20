from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Mapping
from collections.abc import Set
import networkx as nx
import obonet

from .types import GoId, Subontology


_NAMESPACE_TO_SUBONTOLOGY: dict[str, Subontology] = {
    "molecular_function": "MF",
    "biological_process": "BP",
    "cellular_component": "CC",
}

_SUBONTOLOGY_ROOTS: dict[Subontology, GoId] = {
    "MF": "GO:0003674",
    "BP": "GO:0008150",
    "CC": "GO:0005575",
}


@dataclass(frozen=True, slots=True)
class GOTerm:
    """One canonical GO term parsed from an OBO release.

    Attributes
    ----------
    go_id:
        Canonical GO identifier.
    name:
        Human-readable GO term name.
    namespace:
        Raw GO namespace string such as `biological_process`.
    definition:
        Raw GO definition text from the OBO file when present.
    alt_ids:
        Alternative GO identifiers that map to `go_id`.
    parent_ids:
        Canonical parent GO IDs connected by `is_a` or allowed relationships.
    """

    go_id: GoId
    name: str
    namespace: str
    definition: str = ""
    alt_ids: tuple[GoId, ...] = ()
    parent_ids: tuple[GoId, ...] = ()


@dataclass(frozen=True, slots=True)
class GeneOntology:
    """Immutable in-memory GO ontology view.

    Attributes
    ----------
    terms:
        Canonical GO terms keyed by canonical GO ID.
    alt_id_to_term_id:
        Mapping from `alt_id` values to canonical GO IDs.
    release:
        GO release string extracted from the OBO header, typically
        `YYYY-MM-DD`.
    roots:
        Fixed root GO term IDs per subontology.
    """

    terms: dict[GoId, GOTerm]
    alt_id_to_term_id: dict[GoId, GoId]
    release: str | None
    roots: dict[Subontology, GoId]


def read_go_obo(obo_path: str | Path) -> GeneOntology:
    """Parse a GO OBO file into a canonical `GeneOntology` using `obonet`.

    Behavior
    --------
    - keeps only non-obsolete terms
    - resolves `alt_id` mappings to canonical GO IDs
    - keeps `is_a` and `part_of` parent links
    - ignores references to parents that are not present as canonical terms

    Parameters
    ----------
    obo_path:
        Local path to `go-basic.obo`.

    Returns
    -------
    GeneOntology
        Parsed ontology with canonical terms and alternative-ID mapping.
    """

    path = Path(obo_path)
    graph = obonet.read_obo(path)
    release = None
    data_version = graph.graph.get("data-version")
    if isinstance(data_version, str) and data_version.startswith("releases/"):
        release = data_version.split("/", 1)[1]

    canonical_ids = {
        str(node_id)
        for node_id, data in graph.nodes(data=True)
        if not bool(data.get("is_obsolete"))
    }
    terms: dict[GoId, GOTerm] = {}
    alt_id_to_term_id: dict[GoId, GoId] = {}

    for node_id, data in graph.nodes(data=True):
        if str(node_id) not in canonical_ids:
            continue
        go_id = str(node_id)
        raw_alt_ids = data.get("alt_id", [])
        if isinstance(raw_alt_ids, str):
            raw_alt_ids = [raw_alt_ids]
        alt_ids = tuple(str(value) for value in raw_alt_ids)
        parent_ids = tuple(
            str(parent_id)
            for parent_id in graph.successors(node_id)
            if str(parent_id) in canonical_ids
        )
        definition = data.get("def")
        if isinstance(definition, list):
            definition = definition[0] if definition else ""
        parent_ids = tuple(
            parent_id for parent_id in parent_ids if parent_id in canonical_ids
        )
        terms[go_id] = GOTerm(
            go_id=go_id,
            name=str(data.get("name", "")),
            namespace=str(data.get("namespace", "")),
            definition=str(definition or ""),
            alt_ids=alt_ids,
            parent_ids=parent_ids,
        )
        for alt_id in alt_ids:
            alt_id_to_term_id[alt_id] = go_id

    return GeneOntology(
        terms=terms,
        alt_id_to_term_id=alt_id_to_term_id,
        release=release,
        roots=dict(_SUBONTOLOGY_ROOTS),
    )


def canonicalize_go_id(ontology: GeneOntology, term_id: GoId) -> GoId:
    """Return the canonical GO ID for a term or `alt_id`.

    Unknown IDs are returned unchanged so caller code can decide whether to
    ignore or reject them.
    """

    return ontology.alt_id_to_term_id.get(term_id, term_id)


def parents_of(ontology: GeneOntology, term_id: GoId) -> tuple[GoId, ...]:
    """Return direct canonical parents of a GO term.

    Unknown IDs yield an empty tuple.
    """

    canonical_id = canonicalize_go_id(ontology, term_id)
    term = ontology.terms.get(canonical_id)
    if term is None:
        return ()
    return term.parent_ids


def ancestors_of(ontology: GeneOntology, term_id: GoId) -> frozenset[GoId]:
    """Return all recursive canonical ancestors of a GO term."""

    canonical_id = canonicalize_go_id(ontology, term_id)
    if canonical_id not in ontology.terms:
        return frozenset()
    graph = nx.DiGraph(
        (term.go_id, parent_id)
        # (a,b) becomes a --> b.
        # So, the parent here becomes the descendent in the initialzed directed graph
        for term in ontology.terms.values()
        for parent_id in term.parent_ids
    )
    return frozenset(str(node_id) for node_id in nx.descendants(graph, canonical_id))


def propagate_terms(ontology: GeneOntology, term_ids: Set[GoId]) -> frozenset[GoId]:
    """Propagate direct GO terms upward to all canonical ancestors.

    Unknown GO IDs are ignored.
    """

    propagated: set[GoId] = set()
    for term_id in term_ids:
        canonical_id = canonicalize_go_id(ontology, term_id)
        if canonical_id not in ontology.terms:
            continue
        propagated.add(canonical_id)
        propagated.update(ancestors_of(ontology, canonical_id))
    return frozenset(propagated)


def propagate_scores(
    ontology: GeneOntology,
    scores: Mapping[GoId, float],
) -> dict[GoId, float]:
    """Propagate GO scores upward using max child-to-parent semantics.

    Unknown GO IDs are ignored.
    """

    propagated: dict[GoId, float] = {}
    for term_id, score in scores.items():
        canonical_id = canonicalize_go_id(ontology, term_id)
        if canonical_id not in ontology.terms:
            continue
        propagated[canonical_id] = max(float(score), propagated.get(canonical_id, float("-inf")))
        for ancestor_id in ancestors_of(ontology, canonical_id):
            propagated[ancestor_id] = max(float(score), propagated.get(ancestor_id, float("-inf")))
    return propagated


def subontology_terms(ontology: GeneOntology, subontology: Subontology) -> frozenset[GoId]:
    """Return canonical GO terms belonging to one subontology."""

    return frozenset(
        term_id
        for term_id, term in ontology.terms.items()
        if _NAMESPACE_TO_SUBONTOLOGY.get(term.namespace) == subontology
    )


def terms_of_interest(ontology: GeneOntology, subontology: Subontology) -> tuple[GoId, ...]:
    """Return sorted canonical non-root GO terms for one subontology."""

    root_id = ontology.roots[subontology]
    return tuple(sorted(term_id for term_id in subontology_terms(ontology, subontology) if term_id != root_id))


def filter_terms_to_subontology(
    ontology: GeneOntology,
    term_ids: Set[GoId],
    subontology: Subontology,
) -> frozenset[GoId]:
    """Filter a GO-term set down to canonical terms in one subontology."""

    allowed_terms = subontology_terms(ontology, subontology)
    filtered = {
        canonical_id
        for canonical_id in (canonicalize_go_id(ontology, term_id) for term_id in term_ids)
        if canonical_id in allowed_terms
    }
    return frozenset(filtered)
