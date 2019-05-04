#!/usr/bin/env python

"""esgimpl: Concrete implementation for EsgBase.
"""

from esgbase import EsgBase
import augmented_splice_graph as sg

class EsgImpl(EsgBase):
    def load_from_file(self, filename):
        self._splice_graph = sg.load_splice_graph(filename)

    def _unordered_edges(self):
        for edge in self._splice_graph.edges_iter(keys=True, data=True):
            if 'seq' in edge[3]:
                yield edge

    def edges(self):
        return sorted(self._unordered_edges(), key=lambda x: x[0])

    def get_weight(self, e):
        return e[3]['weight']
    
    def get_sequence(self, e, orientation='+'):
        seq = e[3]['seq']
        if orientation == '+':
            return seq
        elif orientation == '-':
            return sg.reverse_complement(seq)
        else:
            raise RuntimeError('Must specify either + or - as orientation')

    def get_orientation(self, e):
        strand = e[3]['strand']
        if strand == sg.Strand.PLUS:
            return '+'
        elif strand == sg.Strand.MINUS:
            return '-'
        else:
            return '*'

    def get_name(self):
        return self._splice_graph.name

    def _unordered_next_edges(self, e, orientation='+'):
        if orientation == '+':
            return self._splice_graph.out_edges_iter([e[1]], keys=True, data=True)
        elif orientation == '-':
            return self._splice_graph.in_edges_iter([e[0]], keys=True, data=True)
        else:
            raise RuntimeError('Must specify either + or - as orientation')

    def get_next_edge(self, e, orientation='+'):
        for edge in sorted(self._unordered_next_edges(e, orientation), key=lambda x: x[3]['weight'], reverse=True):
            if 'seq' in edge[3]:
                yield edge
            else:
                for next_edge in self.get_next_edge(edge, orientation):
                    yield next_edge

    def _to_gene(self, e):
        return sg.find_genes(self._splice_graph, lambda x: e in x.edges_iter(keys=True)).next()
