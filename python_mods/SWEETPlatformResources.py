
from InfoError import *

__all__ = ['SWEETPlatformResources']

class SWEETPlatformResources(InfoError):
	"""
	Information on particular hardware for the target platform (e.g. cluster)

	After each platform sets up the required information, all other variables
	which are still initialized with None are inferred automatically in the
	setup() method.

	In case that insufficient information was provided, an exception is triggered.
	"""

	def __init__(self, dummy_init = False):
		InfoError.__init__(self, "SWEETPlatformResources")

		# Number of physical cores per MPI shared-memory node
		self.num_cores_per_node = None

		# Number of nodes on system
		self.num_nodes = None

		# Total number of physical cores across entire system without hyperthreading
		self.num_cores = None

		# Socket information: Number of physical cores per socket
		self.num_cores_per_socket = None

		# Socket information: Number of sockets per node
		self.num_sockets_per_node = None

		# Maximum wallclock seconds
		self.max_wallclock_seconds = None

		# Number of oversubscribed cores
		# These are the total logically available cores
		# For hyperthreading, these are e.g. the physical ones plus the hyperthreaded ones
		self.num_oversubscribed_cores_per_socket = None


	def setup(self):
		"""
		Setup information after some information was already provided
		"""
		if self.num_cores == None:
			self.num_cores = self.num_cores_per_node*self.num_nodes

		if self.num_sockets_per_node == None:
			if self.num_cores_per_socket == None:
				self.print()
				self.error("Either self.num_sockets_per_node or self.num_cores_per_socket must be specified")
			self.num_sockets_per_node = self.num_cores_per_node // self.num_cores_per_socket

		if self.num_cores_per_socket == None:
			if self.num_sockets_per_node == None:
				self.print()
				self.error("Either self.num_sockets_per_node or self.num_cores_per_socket must be specified")
			self.num_cores_per_socket = self.num_cores_per_node // self.num_sockets_per_node
		if self.num_cores_per_socket*self.num_sockets_per_node != self.num_cores_per_node:
			self.print()
			self.error("Wrong hardware information specified: self.self.num_cores_per_socket*self.num_sockets_per_node != self.num_cores_per_node")

		#
		# Check for all variables being setup
		#
		if self.num_cores_per_node == None:
			self.error("self.num_cores_per_node == None")

		if self.num_nodes == None:
			self.error("self.num_nodes == None")

		if self.num_cores == None:
			self.error("self.num_cores == None")

		if self.num_cores_per_socket == None:
			self.error("self.num_cores_per_socket == None")

		if self.num_sockets_per_node == None:
			self.error("self.num_sockets_per_node == None")


		# Initialize oversubscribed cores per default with number of cores per socket
		if self.num_oversubscribed_cores_per_socket == None:
			self.num_oversubscribed_cores_per_socket = self.num_cores_per_socket



	def print(self):
		self.hline()
		self.info("num_cores_per_node: "+str(self.num_cores_per_node))
		self.info("num_nodes: "+str(self.num_nodes))
		self.info("num_cores: "+str(self.num_cores))
		self.info("num_cores_per_socket: "+str(self.num_cores_per_socket))
		self.info("num_sockets_per_node: "+str(self.num_sockets_per_node))
		self.info("num_oversubscribed_cores_per_socket: "+str(self.num_oversubscribed_cores_per_socket))
		self.info("max_wallclock_seconds: "+str(self.max_wallclock_seconds))



if __name__ == "__main__":
	p = SWEETPlatformResources()
	p.num_cores_per_node = 16
	p.num_cores_per_socket = 4
	p.num_nodes = 12
	p.setup()
	p.print()
	p.info("FIN")
