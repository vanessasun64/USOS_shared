import matplotlib.pyplot as plt
import os
class Dog:
   def __init__(self, name, color, breed):
       self.name = name
       self.color = color
       self.breed = breed
       self.plot_dir = f'{name}_{color}'
       if not os.path.exists(self.plot_dir):
           os.mkdir(self.plot_dir)
   def write_name_into_file(self, fname):
       with open(fname, 'w') as f:
           f.write(self.name)
   def write_name_into_plot(self, fname=None):
       if fname is None:
           fname = os.path.join(self.plot_dir, 'name_plot.png')
       assert fname[-4:] == '.png'
       plt.close()
       plt.xlim(-1, 1); plt.ylim(-1, 1)
       plt.text(0, 0, self.name, va='center', ha='center')
       plt.savefig(fname)
       plt.show()
   def write_breed_into_plot(self, fname=None):
       if fname is None:
           fname = os.path.join(self.plot_dir, 'breed_plot.png')
       assert fname[-4:] == '.png'
       plt.close()
       plt.xlim(-1, 1); plt.ylim(-1, 1)
       plt.text(0, 0, self.breed, va='center', ha='center')
       plt.savefig(fname)
       plt.show()