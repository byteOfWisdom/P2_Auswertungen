from matplotlib import pyplot as plt
import numpy as np


def pd_presets():
	plt.axis('scaled')
	plt.grid()
	plt.xlim([-0.1, 2.5])
	plt.ylim([-0.1, 2.5])
	plt.xlabel("Im(Z)")
	plt.ylabel("Re(Z)")


def pointer(x, y, dx, dy):
	plt.arrow(
		float(x), float(y), 
		float(dx), float(dy), 
		width=0.001, 
		#head_width=0.1, 
		#head_length=0.1, 
		length_includes_head=True,
		color="red"
	)
