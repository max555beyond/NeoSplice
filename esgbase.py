import abc

class EsgBase(object):
	__metaclass__ = abc.ABCMeta

	@abc.abstractmethod
	def load_from_file(self, filename):
		"""Loads a graph from a file.

		Args:
			filename (str): location of the file that holds the graph
		"""
		return

	@abc.abstractmethod
	def edges(self):
		"""Generator function: yields a new edge of the graph
		""" 
		return

	@abc.abstractmethod
	def get_weight(self, e):
		"""Returns the weight of an edge e.
		"""
		return

	@abc.abstractmethod
	def get_sequence(self, e, orientation='+'):
		"""Returns the sequence associated to an edge, in a given orientation.

		Args:
			e: edge of the graph
			orientation ('+' or '-'): the orientation of the sequence. If 
			'+' is specified, the sequence must be returned in transcription
			order. If '-' is specified, the sequence must be returned in reverse
			transcription order.
		"""
		return

	@abc.abstractmethod
	def get_orientation(self, e):
		"""Returns the orientation of an edge ('+' or '-'), or '*' if the
        orientation is ambiguous.

		Args:
			e: edge of the graph
		"""
		return

	@abc.abstractmethod
	def get_name(self):
		"""Returns a name that uniqueliy identifies the graph.
		"""
		return

	@abc.abstractmethod
	def get_next_edge(self, e, orientation='+'):
		"""Generator function. Yields an edge in the graph adjacent to 
		edge e.

		Args:
			e: An edge in the graph
			orientation ('+' or '-'): If '+' is specified, yields edges
			that follow e in the graph in transcription order. If '-' is
			specified, yields edges that precede e in transcription order.
		""" 
		return
