from __future__ import annotations

import copy
import sys
from collections import defaultdict
from dataclasses import dataclass, field
from typing import Callable, Self

COMPLEMENT = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"}


@dataclass
class Segment:
    name: str
    sequence: str

    walks: dict[Walk, int] = field(default_factory=lambda: defaultdict(lambda: 0))
    contexts: dict[tuple[Segment | None, Segment | None], set[Walk]] = field(
        default_factory=lambda: defaultdict(set)
    )

    def __len__(self) -> int:
        return len(self.sequence)

    @property
    def unique_contexts(self) -> int:
        return len([k for k, v in self.contexts.items() if v])

    # @property
    # def unique_contexts(self) -> int:
    #     seen_p, seen_n = set(), set()
    #     cnt = 0
    #     for (p, n), active_walks in self.contexts.items():
    #         if not active_walks:
    #             continue
    #         if p not in seen_p and n not in seen_n:
    #             cnt += 1
    #             seen_p.add(p)
    #             seen_n.add(n)
    #     return cnt

    @property
    def depth(self):
        return sum(v for v in self.walks.values())

    @property
    def duplicated(self):
        return any(v > 1 for v in self.walks.values())

    def __eq__(self, o):
        if not isinstance(o, Segment):
            return False

        return self.name == o.name

    _hash: int | None = None

    def __hash__(self):
        if self._hash is None:
            self._hash = self.name.__hash__()

        return self._hash

    def serialize(self) -> str:
        return "\t".join(["S", self.name, self.sequence])


@dataclass
class OrientedSegment:
    raw_segment: Segment
    reversed: bool = False

    @property
    def sequence(self):
        if self.reversed:
            return "".join(
                COMPLEMENT.get(b, "N") for b in reversed(self.raw_segment.sequence)
            )
        else:
            return self.raw_segment.sequence

    def reverse(self) -> Self:
        s = copy.copy(self)
        s.reversed = not self.reversed
        return s

    def serialize(self):
        return f"{'<' if self.reversed else '>'}{self.raw_segment.name}"


@dataclass
class Link:
    src: Segment
    src_orientation: bool
    dst: Segment
    dst_orientation: bool
    overlap: str

    def serialize(self) -> str:
        return "\t".join([
            "L",
            self.src.name,
            "+" if self.src_orientation else "-",
            self.dst.name,
            "+" if self.dst_orientation else "-",
            self.overlap,
        ])


@dataclass(eq=False)
class Walk:
    sample_id: str
    hap_index: int
    seq_id: str
    start: int
    end: int
    walk: list[OrientedSegment]

    def filter(self, predicate: Callable[[OrientedSegment], bool]):
        new_walk = []
        for i, os in enumerate(self.walk):
            prev = self.walk[i - 1].raw_segment if i > 0 else None
            nxt = self.walk[i + 1].raw_segment if i < len(self.walk) - 1 else None

            if os.reversed:
                os.raw_segment.contexts[(nxt, prev)].discard(self)
            else:
                os.raw_segment.contexts[(prev, nxt)].discard(self)

            if predicate(os):
                new_walk.append(os)
            else:
                os.raw_segment.walks[self] -= 1
                if os.raw_segment.walks[self] <= 0:
                    del os.raw_segment.walks[self]

        for i, os in enumerate(new_walk):
            prev = new_walk[i - 1].raw_segment if i > 0 else None
            nxt = new_walk[i + 1].raw_segment if i < len(new_walk) - 1 else None

            if os.reversed:
                os.raw_segment.contexts[(nxt, prev)].add(self)
            else:
                os.raw_segment.contexts[(prev, nxt)].add(self)

        self.walk = new_walk

    def serialize(self):
        return "\t".join([
            "W",
            self.sample_id,
            str(self.hap_index),
            self.seq_id,
            str(self.start),
            str(self.end),
            "".join(os.serialize() for os in self.walk),
        ])


class Gfa:
    segments: dict[str, Segment] = {}
    links: list[Link] = []
    walks: list[Walk] = []

    def read(self, file_name: str) -> None:
        with open(file_name, "r") as f:
            self.parse(f.read())

    def parse(self, data: str):
        for line in data.splitlines():
            parts = line.split("\t")
            match parts[0]:
                case "#":
                    pass
                case "H":
                    pass
                case "S":
                    s = Segment(parts[1], parts[2])
                    self.segments[s.name] = s
                case "L":
                    src = self.segments[parts[1]]
                    dst = self.segments[parts[3]]

                    self.links.append(
                        Link(src, parts[2] == "+", dst, parts[4] == "+", parts[5])
                    )
                case "W":
                    sample_id = parts[1]
                    hap_index = int(parts[2])
                    seq_id = parts[3]
                    start = int(parts[4])
                    end = int(parts[5])

                    walk: list[OrientedSegment] = []

                    rev: bool = False
                    id: str = ""
                    for c in parts[6]:
                        if c in (">", "<"):
                            if id != "":
                                walk.append(OrientedSegment(self.segments[id], rev))
                            rev = c == "<"
                            id = ""
                        else:
                            id += c
                    if id != "":
                        walk.append(OrientedSegment(self.segments[id], rev))

                    w = Walk(sample_id, hap_index, seq_id, start, end, walk)
                    self.walks.append(w)

                    for os in walk:
                        os.raw_segment.walks[w] += 1

                    for i, os in enumerate(walk):
                        prev = walk[i - 1].raw_segment if i > 0 else None
                        next = walk[i + 1].raw_segment if i < len(walk) - 1 else None
                        if os.reversed:
                            os.raw_segment.contexts[(next, prev)].add(w)
                        else:
                            os.raw_segment.contexts[(prev, next)].add(w)

                case _:
                    raise NotImplementedError(f"{parts[0]} parsing not implemented")

        if not self.walks:
            raise Exception("No walks found, unable to filter.")

    def write(self, file_name: str) -> None:
        with open(file_name, "w") as f:
            f.writelines(s.serialize() + "\n" for s in self.segments.values())
            f.writelines(l.serialize() + "\n" for l in self.links)
            f.writelines(w.serialize() + "\n" for w in self.walks)

    def filter_atoms(
        self,
        min_depth: int = 10,
        max_length: int = 10000,
        max_unique: int = 1,
        remove_dup: bool = False,
    ) -> set[Segment]:
        out = set[Segment]()
        for s in self.segments.values():
            if remove_dup and s.duplicated:
                out.add(s)
            elif len(s) <= max_length:
                if s.depth < min_depth and s.unique_contexts > max_unique:
                    out.add(s)
        return out

    def compute_in_out(
        self,
    ) -> tuple[dict[Segment, set[Segment]], dict[Segment, set[Segment]]]:
        inn: dict[Segment, set[Segment]] = defaultdict(set)
        out: dict[Segment, set[Segment]] = defaultdict(set)
        for w in self.walks:
            for i, os in enumerate(w.walk):
                if os.reversed:
                    if i > 0:
                        out[os.raw_segment].add(w.walk[i - 1].raw_segment)
                    if i < len(w.walk) - 1:
                        inn[os.raw_segment].add(w.walk[i + 1].raw_segment)
                else:
                    if i > 0:
                        inn[os.raw_segment].add(w.walk[i - 1].raw_segment)
                    if i < len(w.walk) - 1:
                        out[os.raw_segment].add(w.walk[i + 1].raw_segment)
        return inn, out

    def find_high_diverse_pairs(
        self, min_in: int = 3, min_out: int = 3, max_span: int = 70000
    ) -> tuple[int, set[Segment]]:
        contexts: set[tuple[Segment, Segment]] = set()

        todo: set[Segment] = set()

        inn, out = self.compute_in_out()

        for w in self.walks:
            n = len(w.walk)
            ps = [0] * (n + 1)
            for i, os in enumerate(w.walk):
                ps[i + 1] = ps[i] + len(os.raw_segment)

            for i, osa in enumerate(w.walk):
                if len(out[osa.raw_segment]) < min_in:
                    continue

                for j in range(i + 1, n):
                    osb = w.walk[j]
                    if len(inn[osb.raw_segment]) < min_out:
                        continue

                    if ps[j + 1] - ps[i] > max_span:
                        break

                    contexts.add((osa.raw_segment, osb.raw_segment))

                    todo.update(w.walk[k].raw_segment for k in range(i + 1, j))

        return len(contexts), todo

    def context_filter(self, segments: set[Segment]) -> dict[Walk, set[Segment]]:
        todo: dict[Walk, set[Segment]] = defaultdict(set)

        for seg in segments:
            if not seg.contexts:
                continue

            best_ctx, _ = max(seg.contexts.items(), key=lambda item: len(item[1]))

            for ctx, walks in seg.contexts.items():
                if ctx == best_ctx:
                    continue
                for w in walks:
                    todo[w].add(seg)

        for w in self.walks:
            w.filter(lambda os: os.raw_segment not in todo[w])

        return todo

    def iterative_filter(self, iterations: int = 1) -> None:
        for it in range(iterations):
            print(f"=== ITERATION {it + 1} ===")

            todo = self.filter_atoms(min_depth=25, remove_dup=True)  # TODO: parameters

            print(f"Globally removing {len(todo)}")

            for w in self.walks:
                w.filter(lambda os: os.raw_segment not in todo)

        contexts, todo_segments = self.find_high_diverse_pairs()  # TODO: args
        print(f"Found {contexts} candidate pairs")

        todo = self.context_filter(todo_segments)

        print(
            f"Context-based removal: atoms removed in {sum(len(v) for v in todo.values())} segments"
        )

        seen_segments = set[Segment]()
        self.walks = [w for w in self.walks if w.walk]

        new_links = {}
        for w in self.walks:
            for os in w.walk:
                seen_segments.add(os.raw_segment)

            for i in range(len(w.walk) - 1):
                src_os = w.walk[i]
                dst_os = w.walk[i + 1]

                src_dir = not src_os.reversed
                dst_dir = not dst_os.reversed

                link_key = (
                    src_os.raw_segment.name,
                    src_dir,
                    dst_os.raw_segment.name,
                    dst_dir,
                )
                if link_key not in new_links:
                    new_links[link_key] = Link(
                        src_os.raw_segment, src_dir, dst_os.raw_segment, dst_dir, "0M"
                    )

        self.links = list(new_links.values())

        self.segments = {k: v for k, v in self.segments.items() if v in seen_segments}
        self.links = [
            l for l in self.links if l.src in seen_segments and l.dst in seen_segments
        ]


if len(sys.argv) < 3:
    print(f"Usage: {sys.argv[0]} in.gfa out.gfa")
    exit(1)

gfa = Gfa()
gfa.read(sys.argv[1])

gfa.iterative_filter(4)

gfa.write(sys.argv[2])
